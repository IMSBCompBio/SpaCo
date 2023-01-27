#library(testthat)
##source("~/SPACO/R/SCA.R")
#source("~/SPACO/R/LoadPreprocessData.R")
#source("~/SPACO/R/AuxiliaryFunctions.R")
# test_that("Alrogithm works independently of order", {
#   expect_equal(SCA_function(data_clean, neighbourindexmatrix,
#                             PC_criterion = "number", PC_value = 100,
#                             SC_criterion = "percent", SC_value = 1,
#                             testStatType = "C"),
#                SCA_function(data_clean[,ncol(data_clean) + 1 - 1:ncol(data_clean)],
#                             neighbourindexmatrix,
#                             PC_criterion = "number", PC_value = 100,
#                             SC_criterion = "percent", SC_value = 1,
#                             testStatType = "C"))
# })
# test_that("First n Spacos same for different SC output values", {
#   expect_equal(SCA_function(data_clean, neighbourindexmatrix,
#                             PC_criterion = "number", PC_value = 100,
#                             SC_criterion = "number", SC_value = 50,
#                             testStatType = "C")[[1]][,1:50],
#                SCA_function(data_clean, neighbourindexmatrix,
#                             PC_criterion = "number", PC_value = 100,
#                             SC_criterion = "number", SC_value = 100,
#                             testStatType = "C")[[1]][,1:50])
# })
# test_that("C computation invariant to scaling", {
#   expect_equal(compute_C(data_clean[1,], neighbourindexmatrix),
#                compute_C(scale(data_clean[1,]), neighbourindexmatrix))
# })
#SCA_Result <- SCA_function(data_clean, neighbourindexmatrix, "percent", 0.8,
                      #   "C", "percent", 1, FALSE)
