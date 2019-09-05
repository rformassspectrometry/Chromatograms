test_that(".valid_chrom_backend_data_storage works", {
    expect_match(.valid_chrom_backend_data_storage(c("a", NA)), "not allowed")
    expect_null(.valid_chrom_backend_data_storage(character()))
    expect_null(.valid_chrom_backend_data_storage("b"))
})
