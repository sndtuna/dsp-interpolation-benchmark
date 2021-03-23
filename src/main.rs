fn main() {
    let mut impulse = vec![0f32; 8];
    let oversample_factor = 8;
    let oversample_factor_recip = (oversample_factor as f32).recip();
    let mut impulse_response = vec![0f32; impulse.len() * oversample_factor];
    impulse[3] = 1f32;
    for response_i in 0..impulse_response.len() {
        let float_index = response_i as f32 * oversample_factor_recip; 
        let int_i = float_index as usize;
        let frac_i = float_index - (int_i as f32);
        impulse_response[response_i] = 
                get_sample_interpolated(&impulse, int_i, frac_i);
    }
}

fn get_sample_interpolated(samples:&[f32], int_i :usize, frac_i :f32) -> f32{
    return 0f32;
}