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

trait Interpolator {
    fn display_name(&self) -> &'static str;
    fn impulse_response_file_name(&self) -> &'static str;
    fn get_width_of_local_context_in_samples(&self) -> usize;
    fn get_sample_interpolated_ref_impl(&self, input: &mut [f32], int_i: isize, 
            frac_i: f32) -> f32;
    fn get_sample_interpolated(&self, input: &mut [f32], int_i: isize, 
            frac_i: f32) -> f32 {
        self.get_sample_interpolated_ref_impl(input, int_i, frac_i)
    }
}

fn main() {
    let warmup_interpolator = Box::new(CPUspeedMeasurement);
    let all_interpolators: Vec<Box<dyn Interpolator>> = vec![
            Box::new(Spline1stDegreeC0),
            Box::new(Spline3rdDegreeC1),
            Box::new(Spline5thDegreeC0),
            Box::new(Spline5thDegreeC1),
            Box::new(Spline5thDegreeC2),
            Box::new(Spline7thDegreeC0),
            Box::new(SincTruncatedApprox),
            Box::new(SincTruncated),
            Box::new(SincHann),
    ];
    let interpolators_with_slower_reference: Vec<Box<dyn Interpolator>> = vec![
            Box::new(Spline1stDegreeC0),
            Box::new(Spline3rdDegreeC1),
            Box::new(Spline5thDegreeC0),
            Box::new(Spline7thDegreeC0),
            Box::new(SincTruncated),
            Box::new(SincHann),
    ];
    let largest_window_size = 
            all_interpolators.iter()
            .map(|x| x.get_width_of_local_context_in_samples())
            .max().unwrap();

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
            6, &*warmup_interpolator, false);
        let time_per_task = now.elapsed().as_secs_f64();
        if time_per_task > 0.025f64 {
            break
        }else{
            testlength_samples *= 2;
        }
    }

    println!("\nUsing test length of {} samples...\n", testlength_samples);
    for _i in 0..3 {
        run_interpolation_and_print_performance(&mut noise, 
            &mut interpolated_noise, oversample_factor, largest_window_size, 
            &*warmup_interpolator, false);
    }
    for interp in &all_interpolators {
        run_interpolation_and_print_performance(&mut noise, 
            &mut interpolated_noise, oversample_factor, largest_window_size,
            &**interp, false);
    }
    println!("\n------unoptimized-reference-implementations------");
    for interp in &interpolators_with_slower_reference {
        run_interpolation_and_print_performance(&mut noise, 
            &mut interpolated_noise, oversample_factor, largest_window_size, 
            &**interp, true);
    }
}

fn run_interpolation_and_print_performance(src: &mut [f32], dest: &mut [f32], 
        oversample_factor: u32, largest_window_size: usize, 
        interp: &dyn Interpolator, use_ref_impl: bool) {
    let now = Instant::now(); 
    let generated_samples = dest.len() as f64;
    resample(src, dest, oversample_factor, largest_window_size, interp, use_ref_impl);
    println!("{:<39}{:>6.1} ns/sample.", interp.display_name().to_owned()+":", 
            now.elapsed().as_secs_f64()*1e9/generated_samples);
}

fn resample(src: &mut [f32], dest: &mut[f32], oversample_factor: u32, 
        largest_window_size: usize, 
        interp: &dyn Interpolator, use_ref_impl: bool) {
    let skip_lower_end = (largest_window_size/2-1) * oversample_factor as usize;
    let skip_higher_end = (largest_window_size/2) * oversample_factor as usize;
    let oversample_factor_recip = (oversample_factor as f64).recip();
    for response_i in skip_lower_end..dest.len()-skip_higher_end {
        let fp_i = response_i as f64 * oversample_factor_recip; 
        let int_i = fp_i as usize;
        let frac_i = fp_i - (int_i as f64);
        dest[response_i] = 
                if use_ref_impl {
                    interp.get_sample_interpolated_ref_impl(src, 
                            int_i as isize, frac_i as f32)
                }else{
                    interp.get_sample_interpolated(src, 
                            int_i as isize, frac_i as f32)
                }
    }
}

struct Spline1stDegreeC0;
impl Interpolator for Spline1stDegreeC0 {
    fn display_name(&self) -> &'static str {
        "1st degree(linear), C0-continuous, 2p"
    } 
    fn impulse_response_file_name(&self) -> &'static str {
        "1st(linear)_C0_IR"
    }
    fn get_width_of_local_context_in_samples(&self) -> usize {
        2
    }
    fn get_sample_interpolated_ref_impl(&self, input: &mut [f32], int_i: isize, 
            frac_i: f32) -> f32 {
        let size = self.get_width_of_local_context_in_samples();
        let y = SlidingWindow::new(input, int_i as usize, size);
        let x = frac_i;
        y[0]*(1.0f32-x) + y[1]*x
    }
    fn get_sample_interpolated(&self, input: &mut [f32], int_i: isize, 
            frac_i: f32) -> f32 {
        let size = self.get_width_of_local_context_in_samples();
        let y = SlidingWindow::new(input, int_i as usize, size);
        let x = frac_i;
        let slope = y[1] - y[0];
        y[0] + slope*x
    }
}

struct Spline3rdDegreeC1;
impl Interpolator for Spline3rdDegreeC1 {
    fn display_name(&self) -> &'static str {
        "3rd degree(cubic), C1-continuous, 4p"
    } 
    fn impulse_response_file_name(&self) -> &'static str {
        "3rd(cubic)_C1_IR"
    }
    fn get_width_of_local_context_in_samples(&self) -> usize {
        4
    }
    fn get_sample_interpolated_ref_impl(&self, input: &mut [f32], int_i: isize, 
            frac_i: f32) -> f32 {
        // references:
        // https://dsp.stackexchange.com/a/18273
        // https://hbfs.wordpress.com/2012/07/03/fast-interpolation-interpolation-part-v/
        // https://en.wikipedia.org/wiki/Cubic_Hermite_spline

        // 4 input samples is the minimum sliding window size needed for cubic splines. 
        let size = self.get_width_of_local_context_in_samples();
        let y = SlidingWindow::new(input, int_i as usize, size);
        // set derivatives at the start/end of the interpolated segment using 
        // central differences (Catmul-Rom).
        let mut dy = [0f32; 2];//only two, because edges have no finite differences.
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
    fn get_sample_interpolated(&self, input: &mut [f32], int_i: isize, 
            frac_i: f32) -> f32 {
        let size = self.get_width_of_local_context_in_samples();
        let y = SlidingWindow::new(input, int_i as usize, size);
        let mut dy = [0f32; 2];

        dy[0] = (y[1] - y[-1])*0.5;
        dy[1] = (y[2] - y[0])*0.5; 
        let y_0m1 = y[0] - y[1];
        let dy_0p1 = dy[0] + dy[1];

        let a: f32 =  2.0*y_0m1          +dy_0p1;
        let b: f32 = -3.0*y_0m1  -dy[0]  -dy_0p1;
        let c: f32 =              dy[0]         ;
        let d: f32 =      y[0]                  ;

        let x = frac_i;
        ((((a*x)+b)*x)+c)*x+d
    }
}

struct CPUspeedMeasurement;
impl Interpolator for CPUspeedMeasurement {
    fn display_name(&self) -> &'static str {
        "(warmup)"
    } 
    fn impulse_response_file_name(&self) -> &'static str {
        "(warmup)"
    }
    fn get_width_of_local_context_in_samples(&self) -> usize {
        4
    }
    fn get_sample_interpolated_ref_impl(&self, input: &mut [f32], int_i: isize, 
            frac_i: f32) -> f32 {
        let size = self.get_width_of_local_context_in_samples();
        let y = SlidingWindow::new(input, int_i as usize, size);
        let mut dy = [0f32; 2];

        dy[0] = (y[1] - y[-1])*0.5;
        dy[1] = (y[2] - y[0])*0.5; 
        let y_0m1 = y[0] - y[1];
        let dy_0p1 = dy[0] + dy[1];

        let a: f32 =  2.0*y_0m1          +dy_0p1;
        let b: f32 = -3.0*y_0m1  -dy[0]  -dy_0p1;
        let c: f32 =              dy[0]         ;
        let d: f32 =      y[0]                  ;

        let x = frac_i;
        ((((a*x)+b)*x)+c)*x+d
    }
}

struct Spline5thDegreeC1;
impl Interpolator for Spline5thDegreeC1 {
    fn display_name(&self) -> &'static str {
        "5th degree, C1-continuous, 6p"
    } 
    fn impulse_response_file_name(&self) -> &'static str {
        "5th_C1_IR"
    }
    fn get_width_of_local_context_in_samples(&self) -> usize {
        6
    }
    fn get_sample_interpolated_ref_impl(&self, input: &mut [f32], int_i: isize, 
            frac_i: f32) -> f32 {
        // references:
        // https://splines.readthedocs.io/en/latest/euclidean/catmull-rom-uniform.html

        // 6 input samples is the minimum sliding window size needed for quintic splines. 
        let size = self.get_width_of_local_context_in_samples();
        let y = SlidingWindow::new(input, int_i as usize, size);
        let x = frac_i;
        // polynomial interpolation of degree N can be made by linear interpolating 
        // between two polynomial interpolations of degree (N-1).
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
}

struct Spline5thDegreeC2;
impl Interpolator for Spline5thDegreeC2 {
    fn display_name(&self) -> &'static str {
        "5th degree, C2-continuous, 6p"
    } 
    fn impulse_response_file_name(&self) -> &'static str {
        "5th_C2_IR"
    }
    fn get_width_of_local_context_in_samples(&self) -> usize {
        6
    }
    fn get_sample_interpolated_ref_impl(&self, input: &mut [f32], int_i: isize, 
            frac_i: f32) -> f32 {
        // references:
        // https://splines.readthedocs.io/en/latest/euclidean/catmull-rom-uniform.html
        // https://en.wikipedia.org/wiki/B-spline#Cardinal_B-spline

        // 6 input samples is the minimum sliding window size needed for quintic splines. 
        let size = self.get_width_of_local_context_in_samples();
        let y = SlidingWindow::new(input, int_i as usize, size);
        let x = frac_i;
        // C2-continuous interpolation of degree N can be made by blending  
        // between three polynomial interpolations of degree (N-2).
        // The blending function can limit the level of smoothness 
        // (continuity of derivatives) of the end result so in this case it 
        // needs to have C1-continuity by itself.
        let mut third_order_lagrange = [0f32; 3];
        third_order_lagrange[0] = 
                 y[-2]        *(x+1.0)*x*(x-1.0)*(-1.0/6.0)
                +y[-1]*(x+2.0)        *x*(x-1.0)*(1.0/2.0)
                +y[ 0]*(x+2.0)*(x+1.0)  *(x-1.0)*(-1.0/2.0)
                +y[ 1]*(x+2.0)*(x+1.0)*x        *(1.0/6.0);

        third_order_lagrange[1] = 
                 y[-1]        *x*(x-1.0)*(x-2.0)*(-1.0/6.0)
                +y[ 0]*(x+1.0)  *(x-1.0)*(x-2.0)*(1.0/2.0)
                +y[ 1]*(x+1.0)*x        *(x-2.0)*(-1.0/2.0)
                +y[ 2]*(x+1.0)*x*(x-1.0)        *(1.0/6.0);

        third_order_lagrange[2] = 
                 y[ 0]  *(x-1.0)*(x-2.0)*(x-3.0)*(-1.0/6.0)
                +y[ 1]*x        *(x-2.0)*(x-3.0)*(1.0/2.0)
                +y[ 2]*x*(x-1.0)        *(x-3.0)*(-1.0/2.0)
                +y[ 3]*x*(x-1.0)*(x-2.0)        *(1.0/6.0);

        // B2-spline blending. 
        // The piecewise definition has x ranging from 0 to 3 so almost all 
        // pieces need to be shifted up from the 0 to 1 range used above. 
        let x01 = x;
        let x12 = x+1.0;
        let x23 = x+2.0;
         third_order_lagrange[2]*0.5*(x01*x01)
        +third_order_lagrange[1]*0.5*(-2.0*x12*x12 + 6.0*x12 -3.0)
        +third_order_lagrange[0]*0.5*(x23-3.0)*(x23-3.0)
    }
}

struct Spline5thDegreeC0;
impl Interpolator for Spline5thDegreeC0 {
    fn display_name(&self) -> &'static str {
        "5th degree, C0-continuous, 6p"
    } 
    fn impulse_response_file_name(&self) -> &'static str {
        "5th_C0_IR"
    }
    fn get_width_of_local_context_in_samples(&self) -> usize {
        6
    }
    fn get_sample_interpolated_ref_impl(&self, input: &mut [f32], int_i: isize, 
            frac_i: f32) -> f32 {
        // 6 input samples is the minimum sliding window size needed for quintic splines. 
        let size = self.get_width_of_local_context_in_samples();
        let y = SlidingWindow::new(input, int_i as usize, size);
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
    fn get_sample_interpolated(&self, input: &mut [f32], int_i: isize, 
            frac_i: f32) -> f32 {
        let size = self.get_width_of_local_context_in_samples();
        let y = SlidingWindow::new(input, int_i as usize, size);
        let x = frac_i;
        // calculations are reordered to allow the compiler to reuse results.
        let fifth_order_lagrange: f32 = 
                            ((x+1.0)*(x*((x-1.0)*((x-2.0)*(x-3.0)))))*(-1.0/120.0)*y[-2] 
                +    (x+2.0)        *(x*((x-1.0)*((x-2.0)*(x-3.0))))   *(1.0/24.0)*y[-1]
                +   ((x+2.0)*(x+1.0))  *((x-1.0)*((x-2.0)*(x-3.0)))   *(-1.0/12.0)*y[ 0]
                +  (((x+2.0)*(x+1.0))*x)        *((x-2.0)*(x-3.0))     *(1.0/12.0)*y[ 1]
                + ((((x+2.0)*(x+1.0))*x)*(x-1.0))        *(x-3.0)     *(-1.0/24.0)*y[ 2]
                +(((((x+2.0)*(x+1.0))*x)*(x-1.0))*(x-2.0))            *(1.0/120.0)*y[ 3];

        fifth_order_lagrange
    }
}

struct Spline7thDegreeC0;
impl Interpolator for Spline7thDegreeC0 {
    fn display_name(&self) -> &'static str {
        "7th degree, C0-continuous, 8p"
    } 
    fn impulse_response_file_name(&self) -> &'static str {
        "7th_C0_IR"
    }
    fn get_width_of_local_context_in_samples(&self) -> usize {
        8
    }
    fn get_sample_interpolated_ref_impl(&self, input: &mut [f32], int_i: isize, 
            frac_i: f32) -> f32 {
        // 8 input samples is the minimum sliding window size needed for quintic splines. 
        let size = self.get_width_of_local_context_in_samples();
        let y = SlidingWindow::new(input, int_i as usize, size);
        let x = frac_i;
        let seventh_order_lagrange: f32 = 
                 y[-3]        *(x+2.0)*(x+1.0)*x*(x-1.0)*(x-2.0)*(x-3.0)*(x-4.0)*(-1.0/5040.0)
                +y[-2]*(x+3.0)        *(x+1.0)*x*(x-1.0)*(x-2.0)*(x-3.0)*(x-4.0)*(1.0/720.0)
                +y[-1]*(x+3.0)*(x+2.0)        *x*(x-1.0)*(x-2.0)*(x-3.0)*(x-4.0)*(-1.0/240.0)
                +y[ 0]*(x+3.0)*(x+2.0)*(x+1.0)  *(x-1.0)*(x-2.0)*(x-3.0)*(x-4.0)*(1.0/144.0)
                +y[ 1]*(x+3.0)*(x+2.0)*(x+1.0)*x        *(x-2.0)*(x-3.0)*(x-4.0)*(-1.0/144.0)
                +y[ 2]*(x+3.0)*(x+2.0)*(x+1.0)*x*(x-1.0)        *(x-3.0)*(x-4.0)*(1.0/240.0)
                +y[ 3]*(x+3.0)*(x+2.0)*(x+1.0)*x*(x-1.0)*(x-2.0)        *(x-4.0)*(-1.0/720.0)
                +y[ 4]*(x+3.0)*(x+2.0)*(x+1.0)*x*(x-1.0)*(x-2.0)*(x-3.0)        *(1.0/5040.0);

        seventh_order_lagrange
    }
    fn get_sample_interpolated(&self, input: &mut [f32], int_i: isize, 
            frac_i: f32) -> f32 {
        let size = self.get_width_of_local_context_in_samples();
        let y = SlidingWindow::new(input, int_i as usize, size);
        let x = frac_i;
        // calculations are reordered to allow the compiler to reuse results.
        let seventh_order_lagrange: f32 = 
                              ((x+2.0)*((x+1.0)*(x*((x-1.0)*((x-2.0)*((x-3.0)*(x-4.0)))))))*(-1.0/5040.0)*y[-3] 
                +      (x+3.0)        *((x+1.0)*(x*((x-1.0)*((x-2.0)*((x-3.0)*(x-4.0))))))   *(1.0/720.0)*y[-2]
                +     ((x+3.0)*(x+2.0))        *(x*((x-1.0)*((x-2.0)*((x-3.0)*(x-4.0)))))   *(-1.0/240.0)*y[-1]
                +    (((x+3.0)*(x+2.0))*(x+1.0))  *((x-1.0)*((x-2.0)*((x-3.0)*(x-4.0))))     *(1.0/144.0)*y[ 0]
                +   ((((x+3.0)*(x+2.0))*(x+1.0))*x)        *((x-2.0)*((x-3.0)*(x-4.0)))     *(-1.0/144.0)*y[ 1]
                +  (((((x+3.0)*(x+2.0))*(x+1.0))*x)*(x-1.0))        *((x-3.0)*(x-4.0))       *(1.0/240.0)*y[ 2]
                + ((((((x+3.0)*(x+2.0))*(x+1.0))*x)*(x-1.0))*(x-2.0))        *(x-4.0)       *(-1.0/720.0)*y[ 3]
                +(((((((x+3.0)*(x+2.0))*(x+1.0))*x)*(x-1.0))*(x-2.0))*(x-3.0))              *(1.0/5040.0)*y[ 4];

        seventh_order_lagrange
    }
}

struct SincTruncated;
impl Interpolator for SincTruncated {
    fn display_name(&self) -> &'static str {
        "truncated sinc, 6p"
    } 
    fn impulse_response_file_name(&self) -> &'static str {
        "truncated_sinc_IR"
    }
    fn get_width_of_local_context_in_samples(&self) -> usize {
        6
    }
    fn get_sample_interpolated_ref_impl(&self, input: &mut [f32], int_i: isize, 
            frac_i: f32) -> f32 {
        let size = self.get_width_of_local_context_in_samples();
        let y = SlidingWindow::new(input, int_i as usize, size);
        let t = frac_i;
        let convolution  = if t==0.0f32 {
            y[0]
        }else{
            let mut sum = 0f32;
            for x_isize in y.x_range() {
                let x = x_isize as f32;
                sum += y[x_isize]*f32::sin((t-x)*consts::PI)/((t-x)*consts::PI);
            }
            sum
        };
        convolution
    }
    fn get_sample_interpolated(&self, input: &mut [f32], int_i: isize, 
            frac_i: f32) -> f32 {
        let size = self.get_width_of_local_context_in_samples();
        let y = SlidingWindow::new(input, int_i as usize, size);
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
}

struct SincHann;
impl Interpolator for SincHann {
    fn display_name(&self) -> &'static str {
        "hann windowed sinc, 6p"
    } 
    fn impulse_response_file_name(&self) -> &'static str {
        "hann_windowed_sinc_IR"
    }
    fn get_width_of_local_context_in_samples(&self) -> usize {
        6
    }
    fn get_sample_interpolated_ref_impl(&self, input: &mut [f32], int_i: isize, 
            frac_i: f32) -> f32 {
        let size = self.get_width_of_local_context_in_samples();
        let hann_window_freq = ((size as f32)*0.5f32).recip();
        let y = SlidingWindow::new(input, int_i as usize, size);
        let t = frac_i;
        let convolution  = if t==0.0f32 {
            y[0]
        }else{
            let mut sum = 0f32;
            for x_isize in y.x_range() {
                let x = x_isize as f32;
                let window = 0.5f32 + 0.5f32*f32::cos((t-x)*consts::PI*hann_window_freq);
                sum += y[x_isize]*window*f32::sin((t-x)*consts::PI)/((t-x)*consts::PI);
            }
            sum
        };
        convolution
    }
    fn get_sample_interpolated(&self, input: &mut [f32], int_i: isize, 
            frac_i: f32) -> f32 {
        let size = self.get_width_of_local_context_in_samples();
        let hann_window_freq = ((size as f32)*0.5f32).recip();
        let y = SlidingWindow::new(input, int_i as usize, size);
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
}

struct SincTruncatedApprox;
impl Interpolator for SincTruncatedApprox {
    fn display_name(&self) -> &'static str {
        "truncated sinc(sin approx.), 6p"
    } 
    fn impulse_response_file_name(&self) -> &'static str {
        "truncated_sinc_sin_approx_IR"
    }
    fn get_width_of_local_context_in_samples(&self) -> usize {
        6
    }
    fn get_sample_interpolated_ref_impl(&self, input: &mut [f32], int_i: isize, 
            frac_i: f32) -> f32 {
        let size = self.get_width_of_local_context_in_samples();
        let y = SlidingWindow::new(input, int_i as usize, size);
        let t = frac_i;
        let convolution  = if t==0.0f32 {
            y[0]
        }else{
            let mut sum = 0f32;
            // linear equations for the sin(x*PI) approximation polynomial:  
            // only half of the sin cycle is approximated.
            // The polynomial needs to pass trough 3 points. The derivative at the 
            // outer two points is also constrained, because the sinc formula is
            // sensitive to that.
            //  ax^4   +bx^3   +cx^2   +dx   +e =   0      (at x=0)
            //  ax^4   +bx^3   +cx^2   +dx   +e =   1      (at x=0.5)
            //  ax^4   +bx^3   +cx^2   +dx   +e =   0      (at x=1)
            // 4ax^3  +3bx^2  +2cx^1   +d       =  PI      (at x=0)
            // 4ax^3  +3bx^2  +2cx^1   +d       = -PI      (at x=1)
            // 
            // solved for a,b,c,d:
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use hound;
    const MAX_DEVIATION_FROM_REFERENCE: f32 = 1e-4; // -80dB

    #[test]
    fn write_impulse_response_files() {
        let mut impulse = vec![0f32; 64];
        impulse[32] = 1.0f32;
        let oversample_factor = 32;
        let mut impulse_response = vec![0f32; impulse.len() * oversample_factor as usize];

        let optimized_implementations: Vec<Box<dyn Interpolator>> = vec![
                Box::new(Spline1stDegreeC0),
                Box::new(Spline3rdDegreeC1),
                Box::new(Spline5thDegreeC1),
                Box::new(Spline5thDegreeC0),
                Box::new(Spline5thDegreeC2),
                Box::new(Spline7thDegreeC0),
                Box::new(SincTruncated),
                Box::new(SincHann),
                Box::new(SincTruncatedApprox),
        ];
        let largest_window_size = 
                optimized_implementations.iter()
                .map(|x| x.get_width_of_local_context_in_samples())
                .max().unwrap();

        for opt_im in optimized_implementations {
            resample(&mut impulse, &mut impulse_response, 
                    oversample_factor, largest_window_size, &*opt_im, false);
            write_to_wav(&impulse_response, oversample_factor, &(
                    "impulse-responses/".to_owned()
                    +opt_im.impulse_response_file_name()
                    +".wav"));
        }
    }

    #[test]
    fn compare_optimized_implementation_with_reference() {
        let mut impulse = vec![0f32; 64];
        impulse[32] = 1.0f32;
        let oversample_factor = 32;
        let mut impulse_response_reference = vec![0f32; impulse.len() * oversample_factor as usize];
        let mut impulse_response = vec![0f32; impulse.len() * oversample_factor as usize];

        let interpolators: Vec<Box<dyn Interpolator>> = vec![
                Box::new(Spline3rdDegreeC1),
                Box::new(Spline5thDegreeC1),
                Box::new(Spline5thDegreeC0),
                Box::new(Spline7thDegreeC0),
                Box::new(SincTruncated),
                Box::new(SincHann),
        ];
        let largest_window_size = 
                interpolators.iter()
                .map(|x| x.get_width_of_local_context_in_samples())
                .max().unwrap();

        for interpolator in interpolators {
            resample(&mut impulse, &mut impulse_response_reference, 
                    oversample_factor, largest_window_size, 
                    &*interpolator, true);
            resample(&mut impulse, &mut impulse_response, 
                    oversample_factor, largest_window_size, 
                    &*interpolator, false);
            compare_arrays(&impulse_response_reference, &impulse_response);
        }
    }

    fn compare_arrays(target: &[f32], actual: &[f32]) {
        for (t, a) in target.iter().zip(actual.iter()) {
            assert!( (t-a).abs()<=MAX_DEVIATION_FROM_REFERENCE );
        }
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