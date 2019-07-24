# PARTITIONR: Hierarchical partitioning of diversity in R

## Version 0.0.1
## Last updated: 24 July 2019

## Getting Help
See our github page at https://github.com/partitionr/PARTITIONR

This version of 'PARTITIONR' is the first publically available R version of the standalone software "PARTITION" developed by T.O. Crist and J.A. Veech (see Crist et al. 2003, DOI: 10.1086/378901 for details). It is highly recommended that you consult 'vignette("PARTITIONR")' to familiarize yourself with the syntax of 'PARTITIONR'. 

### Example

# Install master branch from github
library(devtools)
install_github("partitionr/PARTITIONR@master", build_vignette = TRUE)

# Load library
library(PARTITIONR)

# Read vignette
vignette("PARTITIONR")

# Create fake data
divtest.dat = data.frame(Region = rep(c("Hi","Lo"), each = 4),
                         Site = rep(c("A","B", "C","D"), each = 2),
                         Sample = rep(c("a","b"), times = 4),
                         Sp1 = c(1,2,0,4,5,10,0,1),
                         Sp2 = c(1,0,2,3,1,0,1,0),
                         Sp3 = c(3,0,1,1,1,1,0,1),
                         Sp4 = c(2,2,0,0,5,0,3,4),
                         Sp5 = c(0,0,1,12,3,0,3,5),
                         Sp6 = c(1,0,0,1,4,0,6,0),
                         Sp7 = c(0,0,0,0,7,1,1,3),
                         Sp8 = c(3,5,2,0,0,0,3,1)
)


# Store output as an object
indrich <- partition(data = divtest.dat,
                        levels = c("Sample","Site","Region"),
                        low.level = 1,
                        q = 0,
                        method = "ind",
                        perms = 100)

# Get observed, expected (based on randomized data), and significance tests
summary(indrich)

# Run two-tailed test (either larger or smaller than expected)
summary(indrich, p.value = "two-sided")

# Plot additive beta as a line graph
plot_partition(indrich, beta.type = "add", plot.type = "line")