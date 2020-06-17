#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#define PI 3.14

double* dft(double* signal, int n_dft);
double* get_window(char* win_type, int window_size);
double* windowing(double* signal, double* window, int window_size);
double* masking(double *signal, int mask, int start, int end);
double calc_auto_corr(double* signal, int length, int stride);
double* get_auto_corr(double* signal, int sample_num, int stride);
void durbin_algorithm(double* corrs, double* errors, double* reflect_coef, double **filter_coef, int stride);
double* inverse(double* signal, int length);
double* float_zeros(int size);
double* get_transfer_func(double** filter_coef, int n_dft, int stride);
double* get_log(double* arr, int size);
double* idft(double* spec, int n_dft);
double* lowpass_filter(double* cepstrum, int n_dft, int threshold);
double* highpass_filter(double* cepstrum, int n_dft, int threshold);


double* dft(double* signal, int n_dft) {
	double* spec_magn = float_zeros(n_dft);
	double* spec_real = float_zeros(n_dft);
	double* spec_imag = float_zeros(n_dft);

	for (int k = 0; k < n_dft; k++) {
		spec_real[k] = 0.0;
		spec_imag[k] = 0.0;

		for (int n = 0; n < n_dft; n++) {
			spec_real[k] += signal[n] * cos(2 * PI * k * n / (double)n_dft);
			spec_imag[k] -= signal[n] * sin(2 * PI * k * n / (double)n_dft);
		}
		spec_magn[k] = pow(spec_real[k], 2) + pow(spec_imag[k], 2);
		spec_magn[k] = sqrt(spec_magn[k]);
	}

	free(spec_real);
	free(spec_imag);
	return spec_magn;
}


double* get_window(char* win_type, int window_size) {
	double* window = (double*)malloc(sizeof(double) * window_size);

	if (!strcmp(win_type, "hamming")) {
		for (int n = 0; n < window_size; n++)
			window[n] = 0.54 - 0.46 * cos(2 * PI * n / (window_size - 1));
	}
	else if (!strcmp(win_type, "hanning")) {
		for (int n = 0; n < window_size; n++)
			window[n] = 0.5 - 0.5 * cos(2 * PI * n / (window_size - 1));
	}
	else if (!strcmp(win_type, "rectangular")) {
		for (int n = 0; n < window_size; n++)
			window[n] = 1.0;
	}
	else if (!strcmp(win_type, "blackman")) {
		for (int n = 0; n < window_size; n++)
			window[n] = 0.42 - 0.5 * cos(2 * PI * n / (window_size - 1)) + 0.08 * cos(4 * PI * n / (window_size - 1));
	}
	else if (!strcmp(win_type, "sine")) {
		for (int n = 0; n < window_size; n++)
			window[n] = sin(PI * (n + 0.5) / window_size);
	}
	else {
		printf("%s is not Supported\n", win_type);
	}

	return window;
}


double* windowing(double* signal, double* window, int window_size) {
	for (int n = 0; n < window_size; n++)
		signal[n] *= window[n];

	return signal;
}


double* masking(double *signal, int mask, int start, int end) {
	for (int i = start; i < end; i++) {
		signal[i] = mask;
	}

	return signal;
}


double calc_auto_corr(double* signal, int n_dft, int stride) {
	double corr = 0;

	for (int frame_index = 0; frame_index < n_dft; frame_index++) {
		if (stride <= frame_index)
			corr += signal[frame_index - stride] * signal[frame_index];

		else
			continue;
	}

	return corr;
}


double* get_auto_corr(double* signal, int n_dft, int stride) {
	double* corrs = float_zeros(stride);

	for (int i = 0; i < stride + 1; i++) {
		corrs[i] = calc_auto_corr(signal, n_dft, i);
	}

	return corrs;
}


void durbin_algorithm(double* corrs, double* errors, double* reflect_coef, double **filter_coef, int stride) {
	double sum;
	errors[0] = corrs[0];

	for (int i = 1; i <= stride; i++) {
		sum = 0;
		for (int j = 1; j <= i - 1; j++) {
			sum += filter_coef[i - 1][j] * corrs[i - j];
		}

		reflect_coef[i] = (corrs[i] - sum) / errors[i - 1];
		filter_coef[i][i] = reflect_coef[i];

		for (int j = 1; j <= i - 1; j++) {
			filter_coef[i][j] = filter_coef[i - 1][j] - reflect_coef[i] * filter_coef[i - 1][i - j];
		}
		errors[i] = (1 - pow(reflect_coef[i], 2)) * errors[i - 1];
	}

	return;
}


double* get_transfer_func(double** filter_coef, int n_dft, int stride) {
	double* transfer_func = float_zeros(n_dft);
	transfer_func[0] = 1;

	for (int i = 1; i < n_dft; i++) {
		if (i < stride)
			transfer_func[i] = -filter_coef[stride][i + 1];

		else
			transfer_func[i] = 0;
	}

	return transfer_func;
}


double* float_zeros(int size) {
	double* zeros = (double*)malloc(sizeof(double) * size);

	for (int i = 0; i < size; i++)
		zeros[i] = 0.0;

	return zeros;
}


double* inverse(double* signal, int length) {
	double* inversed = float_zeros(length);

	for (int i = 0; i < length; i++)
		inversed[i] = 1 / signal[i];

	return inversed;
}


double* get_log(double* arr, int size) {
	double *log_arr = float_zeros(size);

	for (int i = 0; i < size; i++)
		log_arr[i] = log10(arr[i]);

	return log_arr;
}


double* idft(double* spec, int n_dft) {
	double *cepstrum = float_zeros(n_dft);
	double *real = float_zeros(n_dft);
	double *imag = float_zeros(n_dft);

	for (int n = 0; n < n_dft; n++) {
		real[n] = imag[n] = 0.0;

		for (int k = 0; k < n_dft; k++) {
			real[n] = real[n] + spec[k] * cos(2 * PI * k * n / (float)n_dft);
			imag[n] = imag[n] + spec[k] * sin(2 * PI * k * n / (float)n_dft);
		}

		cepstrum[n] = real[n] - imag[n];
	}

	free(real);
	free(imag);

	return cepstrum;
}


double* lowpass_filter(double* cepstrum, int n_dft, int threshold) {
	double *filtered = float_zeros(n_dft);

	for (int n = 0; n < n_dft; n++) {
		if (n < threshold || n >= n_dft - threshold)
			filtered[n] = cepstrum[n];
		else
			filtered[n] = 0.0;
	}
	return filtered;
}


double* highpass_filter(double* cepstrum, int n_dft, int threshold) {
	double *filtered = float_zeros(n_dft);

	for (int n = 0; n < n_dft; n++) {
		if (n_dft / 2 - threshold - 1 < n && n < n_dft / 2 + threshold)
			filtered[n] = cepstrum[n];
		else
			filtered[n] = 0.0;
	}
	return filtered;
}