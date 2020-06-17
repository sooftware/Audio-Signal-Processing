import wave


class Pcm2Wav(object):
    def __init__(self):
        self.sample_rate = 16000
        self.channels = 1
        self.bits = 16

    def __call__(self, pcm_path, wav_path):
        pcm = None

        with open(pcm_path, 'rb') as pcmdata:
            pcm = pcmdata.read()

        with wave.open(wav_path, 'wb') as wavdata:
            wavdata.setnchannels(self.channels)
            wavdata.setsampwidth(self.bits // 8)
            wavdata.setframerate(self.sample_rate)
            wavdata.writeframes(pcm)


class Wav2Pcm(object):

    def __call__(self, wav_path, pcm_path):
        sample_Rate, bit_rate, pcmdata = self.convert(wav_path)

        with open(pcm_path, 'wb') as pcm:
            pcm.write(pcmdata)

    def _get_field(self, wav, offset, lent):
        """
        Get values for filed. This is only working for fields with byteorder little
        Args :
          wav : the wave file
          offset : which position to start at.
          lent : length of field
        Return :
          Int of the desired field.
        """
        wav.seek(0)
        wav.seek(offset, 0)
        return int.from_bytes(wav.read(lent), byteorder='little')

    def convert(self, wav_in):
        """
        Get the sample rate, bit rate and PCM raw bytes from a wav.
        Args :
          wav_in : wave file, or string with path to wave file.
        Return :
          sample_rate : int representing the wave file sample rate
          bit_rate : int repesenting the wave file bit rate
          pcm : bytes representing the raw sound.
        """
        if type(wav_in) is str:
            wav_file = open(wav_in, 'rb')
        else:
            wav_file = wav_in
        header_size = self._get_field(wav_file, 16, 4)
        sample_rate = self._get_field(wav_file, 24, 4)
        bit_rate = self._get_field(wav_file, 34, 2)
        wav_file.seek(0)

        if header_size == 16:
            data = wav_file.read()[44:]

        elif header_size == 18:
            data = wav_file.read()[46:]

        else:
            print("WAV format unknown")
            exit(1)

        wav_file.close()
        return sample_rate, bit_rate, data
