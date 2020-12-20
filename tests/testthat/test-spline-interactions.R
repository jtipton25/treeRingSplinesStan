context("spline functions")

test_that("make_Q", {
    # check error messages
    expect_error(make_Q(NA, 1), "df must be a positive integer")
    expect_error(make_Q("aaa", 1), "df must be a positive integer")
    expect_error(make_Q(2.2, 1), "df must be a positive integer")
    expect_error(make_Q(-2, 1), "df must be a positive integer")
    expect_error(make_Q(c(1, 1), 1), "df must be a positive integer")
    expect_error(make_Q(1, 2), "phi must be a number with values between -1 and 1.")
    expect_error(make_Q(1, -2), "phi must be a number with values between -1 and 1.")
    expect_error(make_Q(1, NA), "phi must be a number with values between -1 and 1.")
    expect_error(make_Q(1, "aaa"), "phi must be a number with values between -1 and 1.")
    expect_error(make_Q(1, c(1, 1)), "phi must be a number with values between -1 and 1.")

    expect_identical(
        make_Q(2, 1),
        structure(c(2, -1, -1, 0, -1, 2, 0, -1, -1, 0, 2, -1, 0, -1, -1, 2),
                  .Dim = c(4L, 4L), .Dimnames = list(NULL, NULL))
    )

})


test_that("spline_interactions", {

    set.seed(11)
    x1 <- rnorm(10)
    x2 <- rnorm(10)

    expect_equal(
        spline_interactions(x1, x2, 5),
        structure(c(0.0269309886823068, 0, 0, 0, 0, 0, 0, 0, 9.5744165703055e-06,
                    0.0463063372933386, 0.161723868295918, 0.000318292643574712,
                    0, 0, 0.000192428777776179, 0, 0, 0, 0.0178590562519801, 0.0457105214566745,
                    0.0707050626450576, 0.000675285762823176, 0, 0, 0.0104425338340239,
                    0, 0, 0, 0.0317556788580155, 0.00513003045603798, 0.0021571289342838,
                    0.000219695625257941, 0, 0, 0.146214251775004, 0, 0, 0, 0.00829302526870444,
                    0, 0, 0, 0, 0, 0.366851638608272, 0, 0, 0, 0, 0, 0.0618654599823444,
                    0, 0, 0.185804728260049, 0, 0.0607446088576184, 0, 0, 0.000104837723710557,
                    0.310426793349085, 0.37150962485177, 0.115061748400981, 0, 0.0237528581822273,
                    0.000102653642213684, 0.0838403077796863, 0, 0, 0.195552678466406,
                    0.306432584123885, 0.162422600789653, 0.244113592033071, 0, 0.00049852073692603,
                    0.00557070592242196, 0.01212833126282, 0, 0, 0.34771759322475,
                    0.034390517525991, 0.00495532397027694, 0.0794192491360295, 0,
                    0, 0.0779999003356489, 0, 0, 0, 0.0908067750615285, 0, 0, 0,
                    0, 0, 0.195701793101858, 0, 0.320079763647466, 0, 0, 0, 0.0134116396085355,
                    0, 0, 0.310670573035508, 0, 0.2198949357324, 0, 0, 5.08550872623645e-05,
                    0.119103426058945, 0.0805385299169538, 0.13897084111593, 0, 0.0397154266837374,
                    7.83542135901592e-06, 0.303501157348987, 0, 0, 0.0948594473040509,
                    0.117570942351643, 0.0352111402177164, 0.294838829446968, 0,
                    0.000833540268114879, 0.000425204865877842, 0.0439044496908399,
                    0, 0, 0.168672190889274, 0.0131948290194127, 0.00107425078956591,
                    0.0959220592995566, 0, 0, 0.00595363273929298, 0, 0, 0, 0.0440489005896159,
                    0, 0, 0, 0, 0, 0.014937667837212, 0, 0.594975761928489, 0, 0,
                    0, 0, 0, 0, 0.0357732906740967, 0, 0.101990667876839, 0, 0, 4.4532815735172e-08,
                    0, 0, 0.00799104165067372, 0, 0.0045731769479201, 0, 0.140768524915401,
                    0, 0, 8.30665822228921e-05, 0, 0, 0.0169536958071741, 0, 9.59810194074286e-05,
                    0, 0.020363561952073, 0, 0, 0.000147702973308605, 0, 0, 0.00551566907795945,
                    0, 0, 0, 0, 0, 0, 3.85727697835646e-05, 0, 0, 0, 0, 0, 0, 0,
                    0.0849213026933772, 0, 0, 0, 0, 0, 0, 0, 0, 0.00498608463351133,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00688184314873114, 0, 0.0410096466479079,
                    0, 0, 0, 0, 0, 0, 0, 0.000995526801093188, 0, 0.317711596181992,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0.543109217761054, 0, 0, 0, 0, 0,
                    0, 0, 0, 2.31717306667567e-05, 0.098169539409046, 0, 0),
                  .Dim = c(10L,
                           25L))
    )

    expect_error(spline_interactions(x1, x2, -2), "df must be a positive integer")
    expect_error(spline_interactions(x1, x2, 4.5), "df must be a positive integer")
    x1 <- rnorm(20)
    expect_error(spline_interactions(x1, x2, 5), "x1 and x2 must each be of the same length")
    x1 <- matrix(1:20, 10, 2)
    expect_error(spline_interactions(x1, x2, 5), "x1 must be a numeric vector")
    x1 <- rep("aaa", 10)
    expect_error(spline_interactions(x1, x2, 5), "x1 must be a numeric vector")
    x1 <- rnorm(10)
    x2 <- matrix(1:20, 10, 2)
    expect_error(spline_interactions(x1, x2, 5), "x2 must be a numeric vector")
    x2 <- rep("aaa", 10)
    expect_error(spline_interactions(x1, x2, 5), "x2 must be a numeric vector")

})


test_that("spline_interactions_lattice", {
    set.seed(11)

    expect_equal(
        spline_interactions_lattice(5),
        structure(
            list(
                Var1 = c(1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L,
                         1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L),
                Var2 = c(1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 3L,
                         3L, 3L, 3L, 4L, 4L, 4L, 4L, 4L, 5L, 5L, 5L, 5L, 5L)),
            out.attrs = list(dim = c(5L, 5L), dimnames = list(
                Var1 = c("Var1=1", "Var1=2", "Var1=3", "Var1=4", "Var1=5"),
                Var2 = c("Var2=1", "Var2=2", "Var2=3", "Var2=4", "Var2=5"))),
            class = "data.frame", row.names = c(NA, -25L))
    )

    expect_error(spline_interactions_lattice(-1), "df must be a positive integer")
    expect_error(spline_interactions_lattice(4.5), "df must be a positive integer")
    expect_error(spline_interactions_lattice(NA), "df must be a positive integer")
    expect_error(spline_interactions_lattice("aaa"), "df must be a positive integer")
    expect_error(spline_interactions_lattice(c(3, 4)), "df must be a positive integer")

})

test_that("spline_interactions_penalty", {
    set.seed(11)

    expect_equal(
        spline_interactions_penalty(5),
        structure(c(2, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 3, -1, 0, 0, 0, -1, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 3, -1, 0, 0,
                    0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, -1, 3, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, -1, 2, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 3, -1, 0, 0, 0, -1,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1,
                    4, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, -1, 0, 0, 0, -1, 4, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1, 4, -1, 0, 0, 0, -1,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1,
                    3, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, -1, 0, 0, 0, 0, 3, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1, 4, -1, 0, 0, 0, -1, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1, 4,
                    -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    -1, 0, 0, 0, -1, 4, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1, 3, 0, 0, 0, 0, -1, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 3, -1,
                    0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1,
                    0, 0, 0, -1, 4, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, -1, 0, 0, 0, -1, 4, -1, 0, 0, 0, -1, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1, 4, -1, 0,
                    0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0,
                    0, 0, -1, 3, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, -1, 0, 0, 0, 0, 2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1, 3, -1, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0,
                    -1, 3, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, -1, 0, 0, 0, -1, 3, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1, 2), .Dim = c(25L, 25L
                    ), .Dimnames = list(NULL, NULL))
    )

    expect_error(spline_interactions_penalty(-1), "df must be a positive integer")
    expect_error(spline_interactions_penalty(4.5), "df must be a positive integer")
    expect_error(spline_interactions_penalty(NA), "df must be a positive integer")
    expect_error(spline_interactions_penalty("aaa"), "df must be a positive integer")
    expect_error(spline_interactions_penalty(c(3, 4)), "df must be a positive integer")

    expect_error(spline_interactions_penalty(4, 2), "phi must be a number with values between -1 and 1.")
    expect_error(spline_interactions_penalty(4, -2), "phi must be a number with values between -1 and 1.")
    expect_error(spline_interactions_penalty(4, NA), "phi must be a number with values between -1 and 1.")
    expect_error(spline_interactions_penalty(4, "aaa"), "phi must be a number with values between -1 and 1.")
    expect_error(spline_interactions_penalty(4, c(1, 1)), "phi must be a number with values between -1 and 1.")

})
