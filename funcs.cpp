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




// 1.1
// 模拟VHF信号的生成
std::vector<double> generate_vhf_signal(double frequency, double sample_rate, double duration) {
    std::vector<double> signal;
    double dt = 1.0 / sample_rate;
    for (double t = 0; t < duration; t += dt) {
        double sample = sin(2.0 * M_PI * frequency * t);
        signal.push_back(sample);
    }
    return signal;
}

// 模拟信号传输和接收
std::vector<double> transmit_receive_signal(const std::vector<double>& input_signal, double snr_db) {
    std::vector<double> received_signal;
    double snr_linear = pow(10, snr_db / 10.0);
    std::default_random_engine generator;
    std::normal_distribution<double> noise_distribution(0.0, sqrt(1.0 / (2.0 * snr_linear)));
    
    for (double sample : input_signal) {
        double noise = noise_distribution(generator);
        received_signal.push_back(sample + noise);
    }
    
    return received_signal;
}

int main() {
    double sample_rate = 44100.0;  // 采样率
    double duration = 5.0;        // 信号持续时间（秒）
    double frequency = 144e6;     // VHF频率（Hz）
    double snr_db = 20.0;         // 信噪比（dB）

    std::vector<double> vhf_signal = generate_vhf_signal(frequency, sample_rate, duration);
    std::vector<double> received_signal = transmit_receive_signal(vhf_signal, snr_db);

    // 在这里，您可以将接收到的信号进行进一步处理，如解调等。

    // 输出原始信号和接收到的信号（示例中只输出前10个样本）
    std::cout << "原始VHF信号：" << std::endl;
    for (int i = 0; i < 10; ++i) {
        std::cout << vhf_signal[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "接收到的信号：" << std::endl;
    for (int i = 0; i < 10; ++i) {
        std::cout << received_signal[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}

// 1.2
#include <iostream>
#include <gnuradio/top_block.h>
#include <gnuradio/blocks/file_source.h>
#include <gnuradio/blocks/file_sink.h>
#include <gnuradio/blocks/packed_to_unpacked_bb.h>
#include <gnuradio/blocks/unpacked_to_packed_bb.h>
#include <gnuradio/blocks/const_source_b.h>
#include <gnuradio/blocks/multiply_const_ff.h>
#include <gnuradio/analog/sig_source_c.h>
#include <gnuradio/digital/frequency_modulator_fc.h>
#include <gnuradio/digital/frequency_demodulator_fc.h>
#include <gnuradio/filter/fir_filter_ccf.h>
#include <gnuradio/filter/fir_filter_fff.h>
#include <gnuradio/filter/fir_filter_fff.h>
#include <gnuradio/gr_complex.h>

int main(int argc, char* argv[])
{
    // 创建GNU Radio的顶层流图
    gr::top_block_sptr tb = gr::make_top_block("USRPTransceiver");

    // 创建信号生成器（这里用恒定频率作为示例）
    gr::blocks::const_source_b::sptr source = gr::blocks::const_source_b::make(1.0, 0, 0);

    // 创建调制器
    gr::digital::frequency_modulator_fc::sptr modulator = gr::digital::frequency_modulator_fc::make(1.0);

    // 创建解调器
    gr::digital::frequency_demodulator_fc::sptr demodulator = gr::digital::frequency_demodulator_fc::make(1.0);

    // 创建USRP发送和接收块（需要实际的USRP硬件或模拟器）
    // 在这里，我们只是演示流程，实际上需要连接USRP硬件或模拟器
    gr::blocks::file_sink::sptr file_sink = gr::blocks::file_sink::make(sizeof(float), "received_signal.dat", false);
    gr::blocks::file_source::sptr file_source = gr::blocks::file_source::make(sizeof(float), "transmitted_signal.dat", false, true);

    // 连接信号处理块
    tb->connect(source, 0, modulator, 0);
    tb->connect(modulator, 0, file_sink, 0);

    tb->connect(file_source, 0, demodulator, 0);
    tb->connect(demodulator, 0, file_sink, 0);

    // 运行流图
    tb->run();

    return 0;
}


// 2.1
#include <iostream>
#include <vector>
#include <cmath>
#include <random>

// 模拟信道效应
std::vector<double> simulate_channel(const std::vector<double>& input_signal, double snr_db) {
    // 定义多径效应参数
    const int num_paths = 3;   // 多径通道中的路径数量
    const double path_delays[num_paths] = {0.0, 0.1, 0.2};  // 不同路径的延迟（秒）
    const double path_attenuations[num_paths] = {1.0, 0.8, 0.6};  // 不同路径的衰减

    // 模拟多径效应
    std::vector<double> output_signal(input_signal.size(), 0.0);
    for (int i = 0; i < num_paths; ++i) {
        int delay_samples = static_cast<int>(path_delays[i] * sample_rate); // 延迟样本数
        for (size_t j = delay_samples; j < input_signal.size(); ++j) {
            output_signal[j] += input_signal[j - delay_samples] * path_attenuations[i];
        }
    }

    // 添加噪音
    double snr_linear = pow(10, snr_db / 10.0);
    std::default_random_engine generator;
    std::normal_distribution<double> noise_distribution(0.0, sqrt(1.0 / (2.0 * snr_linear)));
    
    for (size_t i = 0; i < output_signal.size(); ++i) {
        double noise = noise_distribution(generator);
        output_signal[i] += noise;
    }
    
    return output_signal;
}

int main() {
    double sample_rate = 44100.0;  // 采样率
    double duration = 5.0;        // 信号持续时间（秒）
    double snr_db = 20.0;         // 信噪比（dB）

    // 生成原始信号（在这里添加您的信号生成代码）

    // 模拟信道
    std::vector<double> received_signal = simulate_channel(original_signal, snr_db);

    // 在这里，您可以进行信号解调或其他信号处理操作

    // 输出原始信号和接收到的信号（示例中只输出前10个样本）
    std::cout << "原始信号：" << std::endl;
    for (int i = 0; i < 10; ++i) {
        std::cout << original_signal[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "接收到的信号：" << std::endl;
    for (int i = 0; i < 10; ++i) {
        std::cout << received_signal[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}


// 2.2
#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <cmath>

// CommunicationChannel 类
class CommunicationChannel {
public:
    CommunicationChannel(double snr_db) : snr_db_(snr_db) {
        // 初始化随机数生成器
        std::random_device rd;
        generator_ = std::mt19937(rd());
        noise_dist_ = std::normal_distribution<double>(0.0, 1.0);
    }

    // QPSK调制
    std::vector<std::complex<double>> qpsk_modulate(const std::vector<int>& data) {
        std::vector<std::complex<double>> modulated_signal;
        for (int i = 0; i < data.size(); i += 2) {
            int symbol = (data[i] << 1) | data[i + 1];
            double phase = (M_PI / 4) * symbol;
            std::complex<double> modulated_symbol(cos(phase), sin(phase));
            modulated_signal.push_back(modulated_symbol);
        }
        return modulated_signal;
    }

    // QPSK解调
    std::vector<int> qpsk_demodulate(const std::vector<std::complex<double>>& received_signal) {
        std::vector<int> demodulated_data;
        for (const auto& symbol : received_signal) {
            double phase = atan2(imag(symbol), real(symbol));
            int symbol_idx = static_cast<int>((phase + M_PI / 4) / (M_PI / 2)) % 4;
            demodulated_data.push_back((symbol_idx >> 1) & 1);
            demodulated_data.push_back(symbol_idx & 1);
        }
        return demodulated_data;
    }

    // 添加信道效应和噪音
    std::vector<std::complex<double>> channel_with_noise(const std::vector<std::complex<double>>& input_signal) {
        std::vector<std::complex<double>> output_signal;
        double snr_linear = pow(10, snr_db_ / 10.0);
        double noise_variance = 1.0 / (2.0 * snr_linear);

        for (const auto& symbol : input_signal) {
            std::complex<double> noise(noise_dist_(generator_), noise_dist_(generator_));
            output_signal.push_back(symbol + std::sqrt(noise_variance) * noise);
        }

        return output_signal;
    }

private:
    double snr_db_;
    std::mt19937 generator_;
    std::normal_distribution<double> noise_dist_;
};

int main() {
    // 初始化通信信道
    CommunicationChannel channel(20.0); // 设置信噪比为20dB

    // 创建随机QPSK数据
    std::vector<int> qpsk_data;
    for (int i = 0; i < 100; ++i) {
        qpsk_data.push_back(rand() % 2);
        qpsk_data.push_back(rand() % 2);
    }

    // QPSK调制
    std::vector<std::complex<double>> modulated_signal = channel.qpsk_modulate(qpsk_data);

    // 传输信号（添加信道效应和噪音）
    std::vector<std::complex<double>> received_signal = channel.channel_with_noise(modulated_signal);

    // QPSK解调
    std::vector<int> demodulated_data = channel.qpsk_demodulate(received_signal);

    // 输出原始数据和解调后的数据
    std::cout << "原始数据：" << std::endl;
    for (int i = 0; i < qpsk_data.size(); ++i) {
        std::cout << qpsk_data[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "解调后的数据：" << std::endl;
    for (int i = 0; i < demodulated_data.size(); ++i) {
        std::cout << demodulated_data[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}


// 2.3
#include <iostream>
#include <vector>

class CommunicationAlgorithm {
public:
    virtual void process(const std::vector<int>& input_data, std::vector<int>& output_data) = 0;
};

class SimpleAlgorithm : public CommunicationAlgorithm {
public:
    void process(const std::vector<int>& input_data, std::vector<int>& output_data) override {
        // 这里可以是您的通信算法实现
        // 作为示例，我们将输入数据复制到输出数据
        output_data = input_data;
    }
};


class FPGASimulator {
public:
    void loadAlgorithm(CommunicationAlgorithm* algorithm) {
        algorithm_ = algorithm;
    }

    void simulate(const std::vector<int>& input_data, std::vector<int>& output_data) {
        if (algorithm_) {
            algorithm_->process(input_data, output_data);
        } else {
            std::cerr << "No algorithm loaded." << std::endl;
        }
    }

private:
    CommunicationAlgorithm* algorithm_ = nullptr;
};

int main() {
    FPGASimulator simulator;

    // 创建一个简单的通信算法实例
    SimpleAlgorithm simple_algorithm;

    // 加载通信算法到模拟器
    simulator.loadAlgorithm(&simple_algorithm);

    // 模拟通信过程
    std::vector<int> input_data = {1, 0, 1, 0, 1};
    std::vector<int> output_data;
    simulator.simulate(input_data, output_data);

    // 查看输出结果
    std::cout << "Input Data: ";
    for (int bit : input_data) {
        std::cout << bit << " ";
    }
    std::cout << std::endl;

    std::cout << "Output Data: ";
    for (int bit : output_data) {
        std::cout << bit << " ";
    }
    std::cout << std::endl;

    return 0;
}

// 3.1
#include <complex>

class Signal {
public:
    Signal(double amplitude, double phase) : amplitude_(amplitude), phase_(phase) {}

    std::complex<double> getComplexValue() const {
        return amplitude_ * std::exp(std::complex<double>(0, phase_));
    }

private:
    double amplitude_;
    double phase_;
};

class Modulator {
public:
    Modulator() {}

    Signal modulateAM(double data) const {
        // AM调制，使用数据来调制信号的幅度
        return Signal(data, 0.0);
    }

    Signal modulateFM(double data) const {
        // FM调制，使用数据来调制信号的相位
        return Signal(1.0, data);
    }

    Signal modulateQAM(double in_phase_data, double quadrature_data) const {
        // QAM调制，使用数据来调制信号的幅度和相位
        return Signal(std::sqrt(in_phase_data * in_phase_data + quadrature_data * quadrature_data), std::atan2(quadrature_data, in_phase_data));
    }
};

class Demodulator {
public:
    Demodulator() {}

    double demodulateAM(const Signal& received_signal) const {
        // AM解调，提取信号的幅度
        return std::abs(received_signal.getComplexValue());
    }

    double demodulateFM(const Signal& received_signal) const {
        // FM解调，提取信号的相位
        return std::arg(received_signal.getComplexValue());
    }

    void demodulateQAM(const Signal& received_signal, double& in_phase_data, double& quadrature_data) const {
        // QAM解调，提取信号的幅度和相位
        in_phase_data = std::abs(received_signal.getComplexValue());
        quadrature_data = std::arg(received_signal.getComplexValue());
    }
};

int main() {
    Modulator modulator;
    Demodulator demodulator;

    // AM调制
    double am_data = 0.7;
    Signal am_signal = modulator.modulateAM(am_data);
    double demodulated_am_data = demodulator.demodulateAM(am_signal);
    std::cout << "Demodulated AM Data: " << demodulated_am_data << std::endl;

    // FM调制
    double fm_data = 1.2;
    Signal fm_signal = modulator.modulateFM(fm_data);
    double demodulated_fm_data = demodulator.demodulateFM(fm_signal);
    std::cout << "Demodulated FM Data: " << demodulated_fm_data << std::endl;

    // QAM调制
    double qam_in_phase_data = 0.8;
    double qam_quadrature_data = 0.6;
    Signal qam_signal = modulator.modulateQAM(qam_in_phase_data, qam_quadrature_data);
    double demodulated_in_phase_data, demodulated_quadrature_data;
    demodulator.demodulateQAM(qam_signal, demodulated_in_phase_data, demodulated_quadrature_data);
    std::cout << "Demodulated QAM In-Phase Data: " << demodulated_in_phase_data << std::endl;
    std::cout << "Demodulated QAM Quadrature Data: " << demodulated_quadrature_data << std::endl;

    return 0;
}


// 3.2
#include <vector>

class Encoder {
public:
    virtual std::vector<int> encode(const std::vector<int>& input_data) = 0;
};

class HammingEncoder : public Encoder {
public:
    std::vector<int> encode(const std::vector<int>& input_data) override {
        // Hamming编码的实现
        // 这里只是一个示例，实际实现需要更复杂的计算
        std::vector<int> encoded_data;
        // ...
        return encoded_data;
    }
};

class ReedSolomonEncoder : public Encoder {
public:
    std::vector<int> encode(const std::vector<int>& input_data) override {
        // Reed-Solomon编码的实现
        // 这里只是一个示例，实际实现需要更复杂的计算
        std::vector<int> encoded_data;
        // ...
        return encoded_data;
    }
};

#include <vector>

class Decoder {
public:
    virtual std::vector<int> decode(const std::vector<int>& received_data) = 0;
};

class HammingDecoder : public Decoder {
public:
    std::vector<int> decode(const std::vector<int>& received_data) override {
        // Hamming解码的实现
        // 这里只是一个示例，实际实现需要更复杂的计算
        std::vector<int> decoded_data;
        // ...
        return decoded_data;
    }
};

class ReedSolomonDecoder : public Decoder {
public:
    std::vector<int> decode(const std::vector<int>& received_data) override {
        // Reed-Solomon解码的实现
        // 这里只是一个示例，实际实现需要更复杂的计算
        std::vector<int> decoded_data;
        // ...
        return decoded_data;
    }
};

int main() {
    // 创建Hamming编码器和解码器
    HammingEncoder hamming_encoder;
    HammingDecoder hamming_decoder;

    // 创建Reed-Solomon编码器和解码器
    ReedSolomonEncoder rs_encoder;
    ReedSolomonDecoder rs_decoder;

    // 输入数据
    std::vector<int> input_data = {0, 1, 0, 1, 1, 0};

    // 编码数据
    std::vector<int> encoded_hamming_data = hamming_encoder.encode(input_data);
    std::vector<int> encoded_rs_data = rs_encoder.encode(input_data);

    // 模拟信道传输（在实际应用中可能会丢失、损坏数据）

    // 解码数据
    std::vector<int> decoded_hamming_data = hamming_decoder.decode(encoded_hamming_data);
    std::vector<int> decoded_rs_data = rs_decoder.decode(encoded_rs_data);

    // 输出结果
    std::cout << "Input Data: ";
    for (int bit : input_data) {
        std::cout << bit << " ";
    }
    std::cout << std::endl;

    std::cout << "Decoded Hamming Data: ";
    for (int bit : decoded_hamming_data) {
        std::cout << bit << " ";
    }
    std::cout << std::endl;

    std::cout << "Decoded Reed-Solomon Data: ";
    for (int bit : decoded_rs_data) {
        std::cout << bit << " ";
    }
    std::cout << std::endl;

    return 0;
}




// 3.3
#include <iostream>
#include <vector>

// 定义信道效应
class Channel {
public:
    virtual void apply(const std::vector<double>& input_signal, std::vector<double>& output_signal) = 0;
};

// 定义调制器类
class Modulator {
public:
    virtual void modulate(const std::vector<int>& input_data, std::vector<double>& modulated_signal) = 0;
};

// 定义解调器类
class Demodulator {
public:
    virtual void demodulate(const std::vector<double>& received_signal, std::vector<int>& output_data) = 0;
};

// 定义编码器类
class Encoder {
public:
    virtual void encode(const std::vector<int>& input_data, std::vector<int>& encoded_data) = 0;
};

// 定义解码器类
class Decoder {
public:
    virtual void decode(const std::vector<int>& received_data, std::vector<int>& decoded_data) = 0;
};

// 定义通信系统类
class CommunicationSystem {
public:
    CommunicationSystem(Modulator* modulator, Demodulator* demodulator, Encoder* encoder, Decoder* decoder, Channel* channel)
        : modulator_(modulator), demodulator_(demodulator), encoder_(encoder), decoder_(decoder), channel_(channel) {}

    void simulate(const std::vector<int>& input_data, std::vector<int>& received_data) {
        // 调制
        std::vector<double> modulated_signal;
        modulator_->modulate(input_data, modulated_signal);

        // 信道传输
        std::vector<double> channel_output_signal;
        channel_->apply(modulated_signal, channel_output_signal);

        // 解调
        std::vector<int> demodulated_data;
        demodulator_->demodulate(channel_output_signal, demodulated_data);

        // 解码
        std::vector<int> decoded_data;
        decoder_->decode(demodulated_data, decoded_data);

        received_data = decoded_data;
    }

private:
    Modulator* modulator_;
    Demodulator* demodulator_;
    Encoder* encoder_;
    Decoder* decoder_;
    Channel* channel_;
};
