test_that(".valid_chrom_backend_data_storage works", {
    expect_match(.valid_chrom_backend_data_storage(c("a", NA)), "not allowed")
    expect_null(.valid_chrom_backend_data_storage(character()))
    expect_null(.valid_chrom_backend_data_storage("b"))
})

# test_that(".values_match_mz works", {
#     pmz <- c(12.4, 15, 3, 12.4, 3, 1234, 23, 5, 12.4, NA, 3)
#     mz <- c(200, 12.4, 3)
#
#     res <- .values_match_mz(pmz, mz)
#     expect_true(all(pmz[res] %in% mz))
#     expect_false(any(pmz[-res] %in% mz))
#
#     pmz <- rev(pmz)
#     res <- .values_match_mz(pmz, mz)
#     expect_true(all(pmz[res] %in% mz))
#     expect_false(any(pmz[-res] %in% mz))
#
#     res <- .values_match_mz(c(NA, NA), mz)
#     expect_identical(res, integer())
#
#     res <- .values_match_mz(pmz, c(NA, 3))
#     expect_true(all(pmz[res] == 3))
# })



