use std::ops::Index;
use std::ops::IndexMut;
use std::ops::Range;
use std::f32::consts;
use std::time::Instant;
use rand::Rng;
use rand::SeedableRng;
use rand::distributions;

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
        let x_upper_exclusive = self.arr.len() as isize + x_lower_inclusive;
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
    let mut rng = rand::rngs::StdRng::seed_from_u64(0u64);
    let mut testlength_samples: usize = 4096;
    let mut noise: Vec<f32>; 
    let oversample_factor = 32;
    let mut interpolated_noise: Vec<f32>;
    // find a test signal length where the benchmark doesn't finish too quickly if 
    // the fastest interpolation method is used.
    loop { 
        noise = (&mut rng).sample_iter(distributions::Standard)
            .take(testlength_samples).collect();
        interpolated_noise = vec![0f32; noise.len() * oversample_factor as usize];
        let now = Instant::now();
        resample(&mut noise, &mut interpolated_noise, oversample_factor, 
                get_sample_interpolated_linear);
        let time_per_task = now.elapsed().as_secs_f64();
        if time_per_task > 0.025f64 {
            break
        }else{
            testlength_samples *= 2;
        }
    }
    println!("\nUsing test length of {} samples...\n", testlength_samples);

    run_interpolation_and_print_performance(&mut noise,&mut interpolated_noise, 
            oversample_factor, "warmup1", 
            get_sample_interpolated_cubic);

    run_interpolation_and_print_performance(&mut noise,&mut interpolated_noise, 
            oversample_factor, "warmup2", 
            get_sample_interpolated_cubic);

    run_interpolation_and_print_performance(&mut noise,&mut interpolated_noise, 
            oversample_factor, "warmup3", 
            get_sample_interpolated_cubic);

    run_interpolation_and_print_performance(&mut noise,&mut interpolated_noise, 
            oversample_factor, "linear", 
            get_sample_interpolated_linear);

    run_interpolation_and_print_performance(&mut noise,&mut interpolated_noise, 
            oversample_factor, "cubic", 
            get_sample_interpolated_cubic);

    run_interpolation_and_print_performance(&mut noise,&mut interpolated_noise, 
            oversample_factor, "quintic", 
            get_sample_interpolated_quintic);

    run_interpolation_and_print_performance(&mut noise,&mut interpolated_noise, 
            oversample_factor, "quintic pure lagrange", 
            get_sample_interpolated_quintic_pure_lagrange);

    run_interpolation_and_print_performance(&mut noise,&mut interpolated_noise, 
            oversample_factor, "truncated sinc", 
            get_sample_interpolated_truncated_sinc);

    run_interpolation_and_print_performance(&mut noise,&mut interpolated_noise, 
            oversample_factor, "hann windowed sinc", 
            get_sample_interpolated_hann_windowed_sinc);

    run_interpolation_and_print_performance(&mut noise,&mut interpolated_noise, 
            oversample_factor, "truncated sinc(sin approx.)", 
            get_sample_interpolated_truncated_sinc_sin_approx);

}

fn run_interpolation_and_print_performance(src: &mut [f32], dest: &mut [f32], 
        oversample_factor: u32, fn_print_name: &str, 
        interp_func: fn(&mut [f32],isize,f32)->f32 ) {
    let now = Instant::now(); 
    let num_samples = dest.len() as u128;
    resample(src, dest, oversample_factor, interp_func);
    println!("{:<28}{:>6}ns/sample.", fn_print_name.to_owned()+":", 
            now.elapsed().as_nanos()/num_samples);
}

fn resample(src: &mut [f32], dest: &mut[f32], oversample_factor: u32,
            interp_func: fn(&mut [f32],isize,f32)->f32 ) {
    let largest_window_size = 6;
    let skip_lower_end = (largest_window_size/2-1) * oversample_factor as usize;
    let skip_higher_end = (largest_window_size/2) * oversample_factor as usize;
    let oversample_factor_recip = (oversample_factor as f64).recip();
    for response_i in skip_lower_end..dest.len()-skip_higher_end {
        let fp_i = response_i as f64 * oversample_factor_recip; 
        let int_i = fp_i as usize;
        let frac_i = fp_i - (int_i as f64);
        dest[response_i] = 
                interp_func(src, int_i as isize, frac_i as f32);
    }
}

fn get_sample_interpolated_linear(input:&mut [f32], int_i :isize, frac_i :f32) -> f32{
    let y = SlidingWindow::new(input, int_i as usize, 2);
    let x = frac_i;
    y[0]*(1.0-x) + y[1]*x
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

fn get_sample_interpolated_truncated_sinc(input:&mut [f32], int_i :isize, frac_i :f32) -> f32{
    const WIDTH_IN_POINTS: usize = 6;
    // fill sliding window. use zero for indeces outside the input samples.
    let y = SlidingWindow::new(input, int_i as usize, WIDTH_IN_POINTS);
    let t = frac_i;
    let convolution  = if t==0.0f32 {
        y[0]
    }else{
        let mut sum = 0f32;
        let t_pi = t*consts::PI;
        let sin_t_pi = f32::sin(t_pi);
        let mut sin_t_pi_x_pi = -sin_t_pi;
        for x_isize in y.x_range() {
            let x = x_isize as f32;
            let x_pi = x*consts::PI;
            sin_t_pi_x_pi = -sin_t_pi_x_pi; // offsets of PI in sin result in sign flips
            sum += y[x_isize]*sin_t_pi_x_pi/(t_pi-x_pi);
        }
        sum
    };
    convolution
}

fn get_sample_interpolated_hann_windowed_sinc(input:&mut [f32], int_i :isize, frac_i :f32) -> f32{
    const WIDTH_IN_POINTS: usize = 6;
    let hann_window_freq = ((WIDTH_IN_POINTS as f32)*0.5f32).recip();
    // fill sliding window. use zero for indeces outside the input samples.
    let y = SlidingWindow::new(input, int_i as usize, WIDTH_IN_POINTS);
    let t = frac_i;
    let convolution  = if t==0.0f32 {
        y[0]
    }else{
        let mut sum = 0f32;
        let t_pi = t*consts::PI;
        let sin_t_pi = f32::sin(t_pi);
        let mut sin_t_pi_x_pi = -sin_t_pi;
        for x_isize in y.x_range() {
            let x = x_isize as f32;
            let x_pi = x*consts::PI;
            let window = 0.5f32 + 0.5f32*f32::cos((t_pi-x_pi)*hann_window_freq);
            sin_t_pi_x_pi = -sin_t_pi_x_pi; // offsets of PI in sin result in sign flips
            sum += y[x_isize]*window*sin_t_pi_x_pi/(t_pi-x_pi);
        }
        sum
    };
    convolution
}

fn get_sample_interpolated_truncated_sinc_sin_approx(input:&mut [f32], int_i :isize, frac_i :f32) -> f32{
    const WIDTH_IN_POINTS: usize = 6;
    // fill sliding window. use zero for indeces outside the input samples.
    let y = SlidingWindow::new(input, int_i as usize, WIDTH_IN_POINTS);
    let t = frac_i;
    let convolution  = if t==0.0f32 {
        y[0]
    }else{
        let mut sum = 0f32;
        let sin_t_pi = {
            let a =  16.0f32-4.0f32*consts::PI;
            let b = -32.0f32+8.0f32*consts::PI;
            let c =  16.0f32-5.0f32*consts::PI;
            let d =   0.0f32+consts::PI;
            let e =   0.0f32;
            a*(t*t)*(t*t) + b*(t*t)*t + c*(t*t) + d*t + e
        };
        let mut sin_t_pi_x_pi = -sin_t_pi; 
        for x_isize in y.x_range() {
            let x = x_isize as f32;
            sin_t_pi_x_pi = -sin_t_pi_x_pi; // offsets of PI in sin result in sign flips
            sum += y[x_isize]*sin_t_pi_x_pi/(t-x);
        }
        sum*consts::FRAC_1_PI
    };
    convolution
}

#[cfg(test)]
mod tests {
    use super::*;
    use hound;

    #[test]
    fn write_impulse_response_files() {
        let mut impulse = vec![0f32; 64];
        impulse[32] = 1.0f32;
        let oversample_factor = 32;
        let mut impulse_response = vec![0f32; impulse.len() * oversample_factor as usize];

        resample(&mut impulse, &mut impulse_response, oversample_factor,
                get_sample_interpolated_linear);
        write_to_wav(&impulse_response, oversample_factor, 
            "linear_IR.wav");

        resample(&mut impulse, &mut impulse_response, oversample_factor,
                get_sample_interpolated_quintic);
        write_to_wav(&impulse_response, oversample_factor, 
            "quintic_IR.wav");

        resample(&mut impulse, &mut impulse_response, oversample_factor,
                get_sample_interpolated_quintic_pure_lagrange);
        write_to_wav(&impulse_response, oversample_factor, 
            "quintic_pure_lagrange_IR.wav");

        resample(&mut impulse, &mut impulse_response, oversample_factor,
                get_sample_interpolated_cubic);
        write_to_wav(&impulse_response, oversample_factor, 
            "cubic_IR.wav");

        resample(&mut impulse, &mut impulse_response, oversample_factor,
                get_sample_interpolated_truncated_sinc);
        write_to_wav(&impulse_response, oversample_factor, 
            "truncated_sinc_IR.wav");

        resample(&mut impulse, &mut impulse_response, oversample_factor,
            get_sample_interpolated_hann_windowed_sinc);
        write_to_wav(&impulse_response, oversample_factor, 
            "hann_windowed_sinc_IR.wav");

        resample(&mut impulse, &mut impulse_response, oversample_factor,
            get_sample_interpolated_truncated_sinc_sin_approx);
        write_to_wav(&impulse_response, oversample_factor, 
            "truncated_sinc_sin_approx_IR.wav");
    }

    fn write_to_wav(v: &[f32], oversample_factor: u32, filename: &str) {
        let spec = hound::WavSpec {
            channels: 1,
            sample_rate: oversample_factor*48000,
            bits_per_sample: 32,
            sample_format: hound::SampleFormat::Float,
        };
        let mut writer = hound::WavWriter::create(filename, spec).unwrap();
        for s in v.iter() {
            writer.write_sample(*s).unwrap();
        }
        println!("file written.");
    }
}