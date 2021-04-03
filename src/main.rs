use std::ops::Index;
use std::ops::IndexMut;
use std::ops::Range;
use hound;
use std::f32::consts;

/// local context around interpolation segment
struct SlidingWindow<'a> {
    arr: &'a mut [f32],
    zero_index: usize,
}
impl SlidingWindow<'_> {
    fn new<'a>(input: &'a mut [f32], zero_pos: usize, size: usize) -> SlidingWindow<'a> {
        let window_middle_segment_start = (size-1)/2;
        let input_start = zero_pos - window_middle_segment_start;
        let input_end = input_start + size;
        SlidingWindow {
            arr: &mut input[input_start..input_end],
            zero_index: window_middle_segment_start,
        }
    }
    fn x_range(&self) -> Range<isize> {
        let x_lower_inclusive = -(self.zero_index as isize);
        let x_upper_exclusive = self.arr.len() as isize - x_lower_inclusive;
        x_lower_inclusive..x_upper_exclusive
    }
}
// interpolation can only be calculated for the middle segment of the sliding window.
// the formulas assume that the interpolated segment has x values between 0 and 1.
// An Index translation will be used to map x values to window indeces. 
impl Index<isize> for SlidingWindow<'_> {
    type Output = f32;
    fn index(&self, i: isize) -> &f32 {
        &self.arr[(i+(self.zero_index as isize)) as usize]
    }
}
impl IndexMut<isize> for SlidingWindow<'_> {
    fn index_mut(&mut self, i: isize) -> &mut f32 {
        &mut self.arr[(i+(self.zero_index as isize)) as usize]
    }
}

fn main() {
    let mut impulse = vec![0f32; 64];
    impulse[32] = 1.0f32;
    let oversample_factor = 16;
    let mut impulse_response = vec![0f32; impulse.len() * oversample_factor];

    resample(&mut impulse, &mut impulse_response, oversample_factor,
            get_sample_interpolated_quintic);
    write_to_wav(&impulse_response,"quintic_IR.wav");

    resample(&mut impulse, &mut impulse_response, oversample_factor,
            get_sample_interpolated_quintic_pure_lagrange);
    write_to_wav(&impulse_response,"quintic_pure_lagrange_IR.wav");

    resample(&mut impulse, &mut impulse_response, oversample_factor,
            get_sample_interpolated_cubic);
    write_to_wav(&impulse_response,"cubic_IR.wav");

    resample(&mut impulse, &mut impulse_response, oversample_factor,
            get_sample_interpolated_truncated_sinc_6point);
    write_to_wav(&impulse_response,"truncated_sinc_6point_IR.wav");

}

fn resample(src: &mut [f32], dest: &mut[f32], oversample_factor: usize,
            interp_func: fn(&mut [f32],isize,f32)->f32 ) {
    let largest_window_size = 6;
    let skip_lower_end = (largest_window_size/2-1) * oversample_factor;
    let skip_higher_end = (largest_window_size/2) * oversample_factor;
    let oversample_factor_recip = (oversample_factor as f32).recip();
    for response_i in skip_lower_end..dest.len()-skip_higher_end {
        let float_index = response_i as f32 * oversample_factor_recip; 
        let int_i = float_index as usize;
        let frac_i = float_index - (int_i as f32);
        dest[response_i] = 
                interp_func(src, int_i as isize, frac_i);
    }
}

fn write_to_wav(v: &[f32], filename: &str) {
    let spec = hound::WavSpec {
        channels: 1,
        sample_rate: 48000,
        bits_per_sample: 32,
        sample_format: hound::SampleFormat::Float,
    };
    let mut writer = hound::WavWriter::create(filename, spec).unwrap();
    for s in v.iter() {
        writer.write_sample(*s).unwrap();
    }
    println!("file written.");
}

fn get_sample_interpolated_cubic(input:&mut [f32], int_i :isize, frac_i :f32) -> f32{
    // references:
    // https://dsp.stackexchange.com/a/18273
    // https://hbfs.wordpress.com/2012/07/03/fast-interpolation-interpolation-part-v/
    // https://en.wikipedia.org/wiki/Cubic_Hermite_spline

    // 4 input samples is the minimum sliding window size needed for cubic splines. 
    // fill sliding window. use zero for indeces outside the input samples.
    let y = SlidingWindow::new(input, int_i as usize, 4);
    // set derivatives at the start/end of the interpolated segment using 
    // central differences (Catmul-Rom).
    let mut dy = [0f32; 2];//two less, because edges have no finite differences.
    dy[0] = (y[1] - y[-1])*0.5;
    dy[1] = (y[2] - y[0])*0.5; 
    // linear equations that need to be satisfied: 
    //  ax^3   +bx^2  +cx  +d = y[0]      (at x=0)
    //  ax^3   +bx^2  +cx  +d = y[1]      (at x=1)
    // 3ax^2  +2bx    +c      = dy[0]     (at x=0)
    // 3ax^2  +2bx    +c      = dy[1]     (at x=1)
    // 
    // solved for a,b,c,d:
    let a: f32 =  2.0*y[0] -2.0*y[1]     +dy[0]    +dy[1];
    let b: f32 = -3.0*y[0] +3.0*y[1] -2.0*dy[0]    -dy[1];
    let c: f32 =                          dy[0]          ;
    let d: f32 =      y[0]                               ;
    // evaluate cubic at x. (frac_i is already in the same range as x)
    let x = frac_i;
    a*x*x*x + b*x*x + c*x + d
}

fn get_sample_interpolated_quintic(input:&mut [f32], int_i :isize, frac_i :f32) -> f32{
    // references:
    // https://splines.readthedocs.io/en/latest/euclidean/catmull-rom-uniform.html

    // 6 input samples is the minimum sliding window size needed for quintic splines. 
    // fill sliding window. use zero for indeces outside the input samples.
    let y = SlidingWindow::new(input, int_i as usize, 6);
    let x = frac_i;
    let mut fourth_order_lagrange = [0f32; 2];
    fourth_order_lagrange[0] = 
             y[-2]        *(x+1.0)*x*(x-1.0)*(x-2.0)*(1.0/24.0)
            +y[-1]*(x+2.0)        *x*(x-1.0)*(x-2.0)*(-1.0/6.0)
            +y[ 0]*(x+2.0)*(x+1.0)  *(x-1.0)*(x-2.0)*(1.0/4.0)
            +y[ 1]*(x+2.0)*(x+1.0)*x        *(x-2.0)*(-1.0/6.0)
            +y[ 2]*(x+2.0)*(x+1.0)*x*(x-1.0)        *(1.0/24.0);

    fourth_order_lagrange[1] = 
            y[-1]        *x*(x-1.0)*(x-2.0)*(x-3.0)*(1.0/24.0)
           +y[ 0]*(x+1.0)  *(x-1.0)*(x-2.0)*(x-3.0)*(-1.0/6.0)
           +y[ 1]*(x+1.0)*x        *(x-2.0)*(x-3.0)*(1.0/4.0)
           +y[ 2]*(x+1.0)*x*(x-1.0)        *(x-3.0)*(-1.0/6.0)
           +y[ 3]*(x+1.0)*x*(x-1.0)*(x-2.0)        *(1.0/24.0);

     fourth_order_lagrange[0]*(1.0-x)
    +fourth_order_lagrange[1]*x
}

fn get_sample_interpolated_quintic_pure_lagrange(input:&mut [f32], int_i :isize, frac_i :f32) -> f32{
    // 6 input samples is the minimum sliding window size needed for quintic splines. 
    // fill sliding window. use zero for indeces outside the input samples.
    let y = SlidingWindow::new(input, int_i as usize, 6);
    let x = frac_i;
    let fifth_order_lagrange: f32 = 
             y[-2]        *(x+1.0)*x*(x-1.0)*(x-2.0)*(x-3.0)*(-1.0/120.0)
            +y[-1]*(x+2.0)        *x*(x-1.0)*(x-2.0)*(x-3.0)*(1.0/24.0)
            +y[ 0]*(x+2.0)*(x+1.0)  *(x-1.0)*(x-2.0)*(x-3.0)*(-1.0/12.0)
            +y[ 1]*(x+2.0)*(x+1.0)*x        *(x-2.0)*(x-3.0)*(1.0/12.0)
            +y[ 2]*(x+2.0)*(x+1.0)*x*(x-1.0)        *(x-3.0)*(-1.0/24.0)
            +y[ 3]*(x+2.0)*(x+1.0)*x*(x-1.0)*(x-2.0)        *(1.0/120.0);

    fifth_order_lagrange
}

fn get_sample_interpolated_truncated_sinc_6point(input:&mut [f32], int_i :isize, frac_i :f32) -> f32{
    // 6 input samples is the minimum sliding window size needed for quintic splines. 
    // fill sliding window. use zero for indeces outside the input samples.
    let y = SlidingWindow::new(input, int_i as usize, 6);
    let t = frac_i;

    let convolution = 
            if t==0.0f32 {
                y[0]
            }else{ 
                consts::FRAC_1_PI * ( 
                    y[-2]*f32::sin(( 2.0+t)*consts::PI)/( 2.0+t)+
                    y[-1]*f32::sin(( 1.0+t)*consts::PI)/( 1.0+t)+
                    y[ 0]*f32::sin(( 0.0+t)*consts::PI)/( 0.0+t)+
                    y[ 1]*f32::sin((-1.0+t)*consts::PI)/(-1.0+t)+
                    y[ 2]*f32::sin((-2.0+t)*consts::PI)/(-2.0+t)+
                    y[ 3]*f32::sin((-3.0+t)*consts::PI)/(-3.0+t)
                )
            };

    convolution 
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Instant;
    use rand::Rng;
    use rand::SeedableRng;
    use rand::distributions;

    #[test]
    fn benchmark_iterpolators() {
        let rng = rand::rngs::StdRng::seed_from_u64(0u64);
        const TESTLENGTH_SAMPLES: usize = 500000;
        let mut noise: Vec<f32> = rng.sample_iter(distributions::Standard)
                .take(TESTLENGTH_SAMPLES).collect();
        let oversample_factor = 32;
        let mut interpolated_noise = vec![0f32; noise.len() * oversample_factor];

        let now = Instant::now();
        resample(&mut noise, &mut interpolated_noise, oversample_factor,
            get_sample_interpolated_cubic);
        println!("warmup: {} ms.", now.elapsed().as_millis());

        let now = Instant::now();
        resample(&mut noise, &mut interpolated_noise, oversample_factor,
            get_sample_interpolated_cubic);
        println!("warmup2: {} ms.", now.elapsed().as_millis());

        let now = Instant::now();
        resample(&mut noise, &mut interpolated_noise, oversample_factor,
            get_sample_interpolated_cubic);
        println!("warmup3: {} ms.", now.elapsed().as_millis());

        let now = Instant::now();
        resample(&mut noise, &mut interpolated_noise, oversample_factor,
                get_sample_interpolated_cubic);
        println!("cubic: {} ms.", now.elapsed().as_millis());

        let now = Instant::now();
        resample(&mut noise, &mut interpolated_noise, oversample_factor,
                get_sample_interpolated_quintic);
        println!("quintic: {} ms.", now.elapsed().as_millis());

        let now = Instant::now();
        resample(&mut noise, &mut interpolated_noise, oversample_factor,
                get_sample_interpolated_quintic_pure_lagrange);
        println!("quintic_pure_lagrange: {} ms.", now.elapsed().as_millis());

        let now = Instant::now();
        resample(&mut noise, &mut interpolated_noise, oversample_factor,
                get_sample_interpolated_truncated_sinc_6point);
        println!("truncated_sinc_6point: {} ms.", now.elapsed().as_millis());
    }
}