# asymmetric_slepian_wolf

This is a simple code to simulate an asymmetric Slepian-Wolf setting, where the compression/encoding is
done with an LDPC code and the decoding uses side information at the receiver.
For details on the theory, please see the report (pdf file).

The workflow for simulation is as follows.

0. The program expects a parity check matrix 
   in sparse column format, given as two .npy files:
   
   - One file containing the column pointers
   - One file containing the row-indeces
   
   You can either provide those yourself, or use some of the
   codes provided by MacKay in his online archive
   (https://www.inference.org.uk/mackay/codes/data.html).
   
   There you can look for some of the matrices of interest in the alist
   format. Use the instructions on the "convert_alist_to_csc.ipynb" notebook
   to convert it into the required files.
   
1. Adjust the parameters in the "sw_test.cpp" file. 

