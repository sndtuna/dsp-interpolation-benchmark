use std::ops::Index;
use std::ops::IndexMut;
use hound;
fn main() {
    let spec = hound::WavSpec {
        channels: 1,
        sample_rate: 48000,
        bits_per_sample: 32,
        sample_format: hound::SampleFormat::Float,
    };

    let mut impulse = vec![0f32; 8];
    impulse[3] = 1.0f32;
    let oversample_factor = 16;
    let oversample_factor_recip = (oversample_factor as f32).recip();
    let mut impulse_response = vec![0f32; impulse.len() * oversample_factor];
    for response_i in 0..impulse_response.len() {
        let float_index = response_i as f32 * oversample_factor_recip; 
        let int_i = float_index as usize;
        let frac_i = float_index - (int_i as f32);
        impulse_response[response_i] = 
                get_sample_interpolated(&impulse, int_i as isize, frac_i);
    }

    let mut writer = hound::WavWriter::create("spline_IR.wav", spec).unwrap();
    for s in impulse_response.iter() {
        writer.write_sample(*s).unwrap();
    }
    println!("file written.");
}

fn get_sample_interpolated(input:&[f32], int_i :isize, frac_i :f32) -> f32{
    // references:
    // https://dsp.stackexchange.com/a/18273
    // https://hbfs.wordpress.com/2012/07/03/fast-interpolation-interpolation-part-v/
    
    // 4 input samples is the minimum sliding window size needed for cubic splines. 
    struct SlidingWindow {
        arr: [f32; 4],
    }
    // interpolation can only be calculated for the middle segment of the sliding window.
    // the formulas assume that the interpolated segment has x values between 0 and 1.
    // An Index translation will be used to map x values to window indeces. 
    impl Index<isize> for SlidingWindow {
        type Output = f32;
        fn index(&self, i: isize) -> &f32 {
            &self.arr[(i+1) as usize]
        }
    }
    impl IndexMut<isize> for SlidingWindow {
        fn index_mut(&mut self, i: isize) -> &mut f32 {
            &mut self.arr[(i+1) as usize]
        }
    }
    let mut y = SlidingWindow{arr: [0f32; 4]};
    // fill sliding window. use zero for indeces outside the input samples.
    for x in -1..3 as isize {
        if x+int_i >= 0 && x+int_i < input.len() as isize {
            y[x] = input[(x+int_i) as usize];
        }else{
            y[x] = 0.0f32;
        }
    }
    // cubic coefficients
    let a: f32 = -0.5*y[-1] +1.5*y[0] -1.5*y[1] +0.5*y[2];
    let b: f32 =      y[-1] -2.5*y[0] +2.0*y[1] -0.5*y[2];
    let c: f32 = -0.5*y[-1]           +0.5*y[1]          ;
    let d: f32 =                 y[0]                    ;
    // evaluate cubic at x. (frac_i is already in the same range as x)
    let x = frac_i;
    a*x*x*x + b*x*x + c*x + d
}