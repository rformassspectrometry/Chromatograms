test_that(".is_class works", {
    expect_true(.is_class(3L, "integer"))
    expect_true(.is_class(Rle(1:4), "integer"))
})
