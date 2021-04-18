## Overview
This is a small experiment to check the performance of interpolation methods for signal processing applications.

I was curious about methods that go beyond linear interpolation, but aren't restricted to a fixed ratio like FIR-filter up/downsampling. I didn't want to just copy paste magic constants from somewhere so I also wrote some small comments if you are curious how the polynomial interpolators are derived.

You can just look at the results from running the bench mark on the dev machine below or you can run the program yourself.

## Instructions
You need to have Rust installed to run following commands at the root directory of the repro:
* `cargo run` to run the benchmark and output the result in the same terminal.
* `cargo test` to (re-)generate the impulse response files.

## Frequency Spectrum
There is no evaluation of frequency spectrums, but you can use the impulse response files that are included in the repro and open them in a program with a spectrum analyzer like Audacity. To get the most accurate spectrum display with Audacity, set the size equal to the length of the files and use a rectangular window (The files have a length that is a power of two for this purpose).
The impulse response files have a high sampling rate, because they represent a use case of interpolating a 48000Hz sample rate source.

## Results
Output on the dev machine (Ryzen 2700X):
```
Using test length of 262144 samples...

(warmup):                           9.8 ns/sample.
(warmup):                           9.5 ns/sample.
(warmup):                           9.6 ns/sample.
linear:                             5.2 ns/sample.
cubic:                              9.6 ns/sample.
quintic:                           16.4 ns/sample.
quintic pure lagrange:             13.8 ns/sample.
truncated sinc, 6p:                22.8 ns/sample.
hann windowed sinc, 6p:            54.0 ns/sample.
truncated sinc(sin approx.), 6p:   17.0 ns/sample.

------unoptimized-reference-implementations------
cubic:                              9.6 ns/sample.
quintic:                           16.4 ns/sample.
quintic pure lagrange:             13.9 ns/sample.
truncated sinc, 6p:                43.4 ns/sample.
hann windowed sinc, 6p:            71.3 ns/sample.
```

