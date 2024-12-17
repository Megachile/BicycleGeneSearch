Adapted from Stern and Han, Gene Structure-Based Homology Search Identifies Highly Divergent Putative Effector Gene Family
https://academic.oup.com/gbe/article/14/6/evac069/6582783#supplementary-data

Colab files:
Classifier construction
https://colab.research.google.com/drive/1GtSkUDC6US1H_WVIfk3WcJUPrGKI-RY8#scrollTo=Z_Mxc27FsbrS
Classifier application
https://colab.research.google.com/drive/1oJAHehrX4ytqT9OZtLKPWNXQAijsixcx#scrollTo=-OrHiRhcgjKt

train_bicycle_classifier.R corresponds to the classifier construction. Before running it, make sure to load in the necessary files:

# First make sure gdown is installed
conda install -c conda-forge gdown

# Then start R and run:
R

# In R:
system("gdown --id 1gUQxr5dp51igDc-vjvYafEDESPgpdS6Q")  # annotation file
system("gdown --id 1zlLXVQN6kavobWqH29GhUY8AB2VMMv7e")  # annotated transcripts
system("gdown --id 19cyJpv9vC8OuOUr3SgE_wYZDLlKEuti6")  # bicycle genes

then run the classifier:
source("train_bicycle_classifier.R")
main()
q()

Then you can apply the classifier:
Rscript scripts/apply_classifier.R



