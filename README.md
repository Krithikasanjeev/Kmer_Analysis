# Kmer_Analysis
This Git repository comprises various R scripts and Python files aimed at analyzing FASTA sequence data, covering a spectrum of analyses such as Shannon entropy, Super Information, Earth Movers Distance, separation accuracy, Kullback Leibler Divergence, and Jensen Shannon Distance.

Each script serves a distinct purpose:

Entropy_ind_files.R: Computes Shannon entropy based on the occurrence frequency of overlapping k-mers across all files, disregarding their profile affiliation.

Entropy_Combined_files.R: Determines Shannon entropy for each profile by tallying the occurrences of overlapping k-mers within that specific profile.

kl_dist.R: Computes Kullback Leibler divergence and Jensen Shannon Distance by evaluating the probability distribution of k-mers across 
different datasets.

Super_information.R: Calculates Super Information by initially computing the Shannon entropy of non-overlapping k-mers, followed by determining the nucleotide occurrences within each k-mer and subsequently computing their entropy, providing a higher-level entropy measure.

Accuracy_of_separation.R: Evaluates the Accuracy of Separation by computing the overlap integral between probability distributions.

Earth_movers_distance: Determines the Earthmovers Distance by assessing the cumulative difference between distributions.

Logistic_Regression.py: Utilizes Logistic Regression models in Python to train datasets for profiling analysis.

Random_forest.py: Utilizes Random Forest models in Python to train datasets for profiling analysis.
These scripts collectively offer a comprehensive toolkit for analyzing sequence data and classifying profiles using both statistical and machine learning techniques.
