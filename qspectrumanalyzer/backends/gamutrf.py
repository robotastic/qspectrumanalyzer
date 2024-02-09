import math, shlex
import zmq
import pandas as pd
import datetime
import json
import logging

from Qt import QtCore

from qspectrumanalyzer import subprocess
from qspectrumanalyzer.backends import BaseInfo, BasePowerThread


class Info(BaseInfo):
    """rtl_power_fftw device metadata"""
    pass

class PowerThread(BasePowerThread):
    """Thread which runs rtl_power_fftw process"""
    def setup(self, start_freq, stop_freq, bin_size, interval=10.0, gain=-1, ppm=0, crop=0,
              single_shot=False, device=0, sample_rate=2560000, bandwidth=0, lnb_lo=0):

        zmq_context = zmq.Context()
        self.socket = zmq_context.socket(zmq.SUB)
        self.socket.connect("tcp://localhost:10000")
        self.socket.setsockopt_string(zmq.SUBSCRIBE, "")


        """Setup rtl_power_fftw params"""
        crop = crop * 100
        overlap = crop * 2
        freq_range = stop_freq * 1e6 - start_freq * 1e6
        min_overhang = sample_rate * overlap * 0.01
        hops = math.ceil((freq_range - min_overhang) / (sample_rate - min_overhang))
        overhang = (hops * sample_rate - freq_range) / (hops - 1) if hops > 1 else 0
        if bin_size > 2800:
            bin_size = 2800
        bins = math.ceil(sample_rate / (bin_size * 1e3))
        crop_freq = sample_rate * crop * 0.01

        self.params = {
            "start_freq": start_freq,
            "stop_freq": stop_freq,
            "freq_range": freq_range,
            "device": device,
            "sample_rate": int(sample_rate),
            "bin_size": bin_size,
            "bins": bins,
            "interval": interval,
            "hops": hops,
            "time": interval / hops,
            "gain": int(gain * 10),
            "ppm": ppm,
            "crop": crop,
            "overlap": overlap,
            "min_overhang": min_overhang,
            "overhang": overhang,
            "single_shot": single_shot
        }
        self.lnb_lo = lnb_lo
        self.freqs = [self.get_hop_freq(hop) for hop in range(hops)]
        self.freqs_crop = [(f[0] + crop_freq, f[1] - crop_freq) for f in self.freqs]
        self.databuffer = {"timestamp": [], "x": [], "y": []}
        self.databuffer_hop = {"timestamp": [], "x": [], "y": []}
        self.hop = 0
        self.prev_line = ""
        self.txt_buf = ""
        self.scan_configs = {}
        self.fftbuffer = None

    def get_hop_freq(self, hop):
        """Get start and stop frequency for particular hop"""
        start_freq = self.params["start_freq"] * 1e6 + (self.params["sample_rate"] - self.params["overhang"]) * hop
        stop_freq = start_freq + self.params["sample_rate"] - (self.params["sample_rate"] / self.params["bins"])
        return (start_freq, stop_freq)


    def process_start(self):
        """Start rtl_power_fftw process"""



    def parse_output(self, line):
        """Parse one line of output from rtl_power_fftw"""
        line = line.strip()

        # One empty line => new hop
        if not line and self.prev_line:
            self.hop += 1
            self.databuffer["x"].extend(self.databuffer_hop["x"])
            self.databuffer["y"].extend(self.databuffer_hop["y"])
            self.databuffer_hop = {"timestamp": [], "x": [], "y": []}

        # Two empty lines => new set
        elif not line and not self.prev_line:
            self.hop = 0
            self.data_storage.update(self.databuffer)
            self.databuffer = {"timestamp": [], "x": [], "y": []}

        # Get timestamp for new hop and set
        elif line.startswith("# Acquisition start:"):
            timestamp = line.split(":", 1)[1].strip()
            if not self.databuffer_hop["timestamp"]:
                self.databuffer_hop["timestamp"] = timestamp
            if not self.databuffer["timestamp"]:
                self.databuffer["timestamp"] = timestamp

        # Skip other comments
        elif line.startswith("#"):
            pass

        # Parse frequency and power
        elif line[0].isdigit():
            freq, power = line.split()
            freq, power = float(freq) + self.lnb_lo, float(power)
            start_freq, stop_freq = self.freqs_crop[self.hop]

            # Apply cropping
            if freq >= start_freq and freq <= stop_freq:
                # Skip overlapping frequencies
                if not self.databuffer["x"] or freq > self.databuffer["x"][-1]:
                    #print("  {:.3f} MHz".format(freq / 1e6))
                    self.databuffer_hop["x"].append(freq)
                    self.databuffer_hop["y"].append(power)
                else:
                    #print("  Overlapping {:.3f} MHz".format(freq / 1e6))
                    pass
            else:
                #print("  Cropping {:.3f} MHz".format(freq / 1e6))
                pass

        self.prev_line = line

    def read_new_frame_df(self, df, discard_time):
        frame_df = None
        scan_config = None
        if discard_time:
            df = df[(time.time() - df.ts).abs() < discard_time]
        if df.size:
            lastfreq = df["freq"].iat[-1]
            #print("last frequency read %f MHz" % (lastfreq / 1e6))
            if self.fftbuffer is None:
                self.fftbuffer = df
            else:
                self.fftbuffer = pd.concat([self.fftbuffer, df])
            if self.fftbuffer["sweep_start"].nunique() > 1:
                min_sweep_start = self.fftbuffer["sweep_start"].min()
                max_sweep_start = self.fftbuffer["sweep_start"].max()
                frame_df = self.fftbuffer[
                    self.fftbuffer["sweep_start"] != max_sweep_start
                ].copy()
                frame_df["tune_count"] = (
                    frame_df["tune_count"].max() - frame_df["tune_count"].min()
                )
                self.fftbuffer = self.fftbuffer[
                    self.fftbuffer["sweep_start"] == max_sweep_start
                ]
                scan_config = self.scan_configs[min_sweep_start]
                del self.scan_configs[min_sweep_start]
        #print("last frequency read %f MHz" % (lastfreq / 1e6))
        return (scan_config, frame_df)

    def lines_to_df(self, lines):
        try:
            records = []
            for line in lines:
                line = line.strip()
                json_record = json.loads(line)
                ts = float(json_record["ts"])
                sweep_start = float(json_record["sweep_start"])
                total_tune_count = int(json_record["total_tune_count"])
                buckets = json_record["buckets"]
                scan_config = json_record["config"]
                self.scan_configs[sweep_start] = scan_config
                records.extend(
                    [
                        {
                            "ts": ts,
                            "freq": float(freq),
                            "db": float(db),
                            "sweep_start": sweep_start,
                            "tune_count": total_tune_count,
                        }
                        for freq, db in buckets.items()
                    ]
                )
            return pd.DataFrame(records)
        except ValueError as err:
            logging.error(str(err))
            return None

    def run(self):
        """hackrf_sweep thread main loop"""
        self.process_start()
        self.alive = True
        self.powerThreadStarted.emit()

        while self.alive:
            self.txt_buf += self.socket.recv().decode("utf-8") #(flags=zmq.NOBLOCK)
            lines = self.txt_buf.splitlines()
            #print(lines)
            if len(lines) > 1:
                if self.txt_buf.endswith("\n"):
                    self.txt_buf = ""
                elif lines:
                    last_line = lines[-1]
                    self.txt_buf = last_line
                    lines = lines[:-1]
                df = self.lines_to_df(lines)
                if df is not None:
                    scan_config, frame_df = self.read_new_frame_df(df, False)   
                    
                    if frame_df is not None:
                        self.databuffer["x"] = frame_df['freq'].to_numpy()
                        self.databuffer["y"] = frame_df['db'].to_numpy()
                        self.data_storage.update(self.databuffer)    



        self.process_stop()
        self.alive = False
        self.powerThreadStopped.emit()