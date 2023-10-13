#include <iostream>
#include <vector>
#include <cmath>
#include <complex>

// ################## hanning window ##################
std::vector<double> hanningWindow(int windowLength) {
    std::vector<double> window;
    for (int n = 0; n < windowLength; ++n) {
        double val = 0.5 - 0.5 * cos((2 * M_PI * n) / (windowLength - 1));
        window.push_back(val);
    }
    return window;
}

int main() {
    int windowLength = 64;
    std::vector<double> hanningWindowValues = hanningWindow(windowLength);

    // Print the Hanning window values
    for (double value : hanningWindowValues) {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    return 0;
}


// ################## (inner) convolution/cross correlation ##################
std::vector<double> innerConvolution(const std::vector<double>& input_signal, const std::vector<double>& kernel) {
    int input_length = input_signal.size();
    int kernel_length = kernel.size();
    int output_length = input_length - kernel_length + 1;

    std::vector<double> output_signal(output_length, 0.0);

    for (int i = 0; i < output_length; ++i) {
        for (int j = 0; j < kernel_length; ++j) {
            output_signal[i] += input_signal[i + j] * kernel[j];
        }
    }

    return output_signal;
}


############ low-pass filter ####################
// Function to apply a low-pass filter to an input signal
std::vector<double> lowPassFilter(const std::vector<double>& input, double cutoffFrequency, double samplingRate) {
    int numSamples = input.size();
    std::vector<double> output(numSamples, 0.0);

    // Calculate the filter's time constant
    double tau = 1.0 / (2.0 * 3.14159265359 * cutoffFrequency);

    // Calculate the filter coefficients
    double alpha = 1.0 / (1.0 + tau * samplingRate);
    for (int i = 1; i < numSamples; ++i) {
        output[i] = alpha * input[i] + (1.0 - alpha) * output[i - 1];
    }

    return output;
}


// ############################## Define the FFT function ####################
void fft(std::vector<std::complex<double>>& x) {
    int N = x.size();
    if (N <= 1) {
        return;
    }

    // Divide the even and odd-indexed elements
    std::vector<std::complex<double>> even(N / 2), odd(N / 2);
    for (int i = 0; i < N / 2; ++i) {
        even[i] = x[i * 2];
        odd[i] = x[i * 2 + 1];
    }

    // Recursively compute FFT for even and odd parts
    fft(even);
    fft(odd);

    // Combine the results
    for (int k = 0; k < N / 2; ++k) {
        std::complex<double> t = std::polar(1.0, -2 * M_PI * k / N) * odd[k];
        x[k] = even[k] + t;
        x[k + N / 2] = even[k] - t;
    }
}



// ###############   IFFT ##################################
void ifft(std::vector<std::complex<double>>& x) {
    int N = x.size();
    if (N <= 1) {
        return;
    }

    // Conjugate the input
    for (int i = 0; i < N; ++i) {
        x[i] = std::conj(x[i]);
    }

    // Compute FFT using the same FFT function
    fft(x);

    // Conjugate the output and scale
    for (int i = 0; i < N; ++i) {
        x[i] = std::conj(x[i]) / static_cast<double>(N);
    }
}



// ######################### Quantization function ###################
std::vector<int> quantize(const std::vector<double>& input, int numLevels) {
    std::vector<int> quantizedValues;
    double minValue = *min_element(input.begin(), input.end());
    double maxValue = *max_element(input.begin(), input.end());

    double stepSize = (maxValue - minValue) / (numLevels - 1);

    for (double value : input) {
        // Quantize the value to the nearest level
        int quantizedValue = static_cast<int>(round((value - minValue) / stepSize));
        quantizedValues.push_back(quantizedValue);
    }

    return quantizedValues;
}


 
// ######################Amplitude Modulation (AM) function ###########################
std::vector<double> amplitudeModulation(const std::vector<double>& message, double carrierFrequency, double modulationIndex, double samplingRate, double duration) {
    int numSamples = static_cast<int>(duration * samplingRate);
    std::vector<double> modulatedSignal(numSamples);

    for (int i = 0; i < numSamples; ++i) {
        double time = static_cast<double>(i) / samplingRate;
        double messageSignal = message[i % message.size()]; // Repeating the message signal if it's shorter

        // Calculate the carrier signal
        double carrierSignal = sin(2 * PI * carrierFrequency * time);

        // Perform AM modulation
        modulatedSignal[i] = (1 + modulationIndex * messageSignal) * carrierSignal;
    }

    return modulatedSignal;
}