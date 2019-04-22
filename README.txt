------------------
Usage Instructions
------------------

Linux instructions:

    To run this program, unzip our project. Change directories
    to our project directory (the directory containing train_and_test.py).

    To train our model, run

        python train_and_test.py <pssm_directory> <tm_align_directory> <fasta_directory>

    where pssm_directory is the directory path containing the pssm files,
    tm_align_directory is the directory containing the tmalign files, and
    fasta_directory is the directory containing the fasta files.

    To predict a new tm-score, run

        python predict_tm_score.py <pssm_directory> <pssm_file_name_1> <pssm_file_name_2> <fasta_directory>

    where pssm_directory is the directory path containing the pssm files,
    pssm_file_name_1 and pssm_file_name_2 are the file names of pssm files within
    pssm_directory (Example: pssm_file_name_1 = '1a3a.pssm'), and
    fasta_directory is the directory containing the fasta files.


Notes:
    Runs on python2 and python3
    Tested on Ubuntu 18.04. Not tested on Windows or Mac.


Examples:

    # Running training & testing
    python3.6 train_and_test.py 5970_6970_SP_19_PROJECT_5/5970_6970_SP_19_PROJECT_5/pssm 5970_6970_SP_19_PROJECT_5/5970_6970_SP_19_PROJECT_5/tmalign 5970_6970_SP_19_PROJECT_5/5970_6970_SP_19_PROJECT_5/fasta

    # Running prediction
    python3.6 predict_tm_score.py 5970_6970_SP_19_PROJECT_5/5970_6970_SP_19_PROJECT_5/pssm 1a3a.pssm 1a6m.pssm 5970_6970_SP_19_PROJECT_5/5970_6970_SP_19_PROJECT_5/fasta
