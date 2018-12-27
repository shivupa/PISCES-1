
void DVR::FFTSetup()
{


    int nthreads = 1;    
    nthreads = omp_get_max_threads();


    int nrofpts = ;
    KE_diag = new double[nrofpts];        
    phi_x       = new Complex[nrofpts];     
    phi_k       = new Complex[nrofpts]; 
    KE_phi_k    = new Complex[nrofpts];
    KE_phi_x    = new Complex[nrofpts];

    printf("FFT Using %u inti threads.\n",fftw_init_threads() );
    fftw_plan_with_nthreads(nthreads);
    printf("FFT Using %u OpenMP threads.\n",omp_get_max_threads() );

    //Define the forward and backward FFT plans
     plan_forward   = fftw_plan_dft(3, n_1dbas,(fftw_complex*)phi_x,(fftw_complex*)phi_k, FFTW_FORWARD,  FFTW_ESTIMATE);
     plan_backward  = fftw_plan_dft(3, n_1dbas,(fftw_complex*)KE_phi_k,(fftw_complex*)KE_phi_x, FFTW_BACKWARD, FFTW_ESTIMATE);

}

