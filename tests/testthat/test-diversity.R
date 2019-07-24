context("Testing diversity measures")
library(PARTITIONR)


divtest.dat = data.frame(Site = rep(c("A","B"), each = 4),
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

partytest0 <- partition(sp = divtest.dat,
                       h.level = c("Sample","Site"),
                       low.level = 1,
                       q = 0,
                       hyp.test = "NONE")
partytest1 <- partition(sp = divtest.dat,
                        h.level = c("Sample","Site"),
                        low.level = 1,
                        q = 1,
                        hyp.test = "NONE")

expect_error(partition(sp = divtest.dat,
                       h.level = c("Site","Sample"),
                       low.level = 1,
                       q = sjlehbrhe,
                       hyp.test = "NONE"),
             "object 'sjlehbrhe' not found")


expect_error(partition(sp = divtest.dat,
                       h.level = c("Site","Sample"),
                       low.level = 1,
                       q = 0,
                       hyp.test = "Test"),
             "hyp.test must be set to \"NONE\", \"INDIVIDUAL\", or \"SAMPLE\"")

expect_equal(as.numeric(partytest0[[1]][1]),8)
expect_equal(as.numeric(partytest1[[1]][1]),7.4170256)

expect_is(partytest0, "partition")
