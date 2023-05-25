# Epoch Based Spectrum Estimation for Speech 
## Description
This repository allows you to regenerate the results from the paper [Gudnason, J., Fang, G., Brookes. M, "Epoch Based Spectrum Estimation for Speech" in Proc. Interspeech, Dublin, Ireland, 2023](https://lvl.ru.is/the-team/publications/). The code produces spectrograms that are in synchrony with epochs in the speech signal derived from glottal glosure instants.  A copy-synthesis scheme using Mel-filter energies is set up for both epoch based and fixed frame spectrograms where SNR and PESQ scores can be compared.

## Setup
The code is developed in Matlab Version: 9.14.0.2206163 (R2023a) and requires [Voicebox](http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html) to be installed.  It also relies on the TIMIT dataset.  In order to run the code, create your own m-file called **projParam.m**.  You can do this for example by copying the content of **projParamDemo.m** and modify it accordingly, for example by entering the root directory of TIMIT.

## License
See the license file for detail.
