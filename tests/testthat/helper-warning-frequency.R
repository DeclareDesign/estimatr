# `check_se_type()` signals its two default-substitution warnings with
# `rlang::warn(.frequency = "once")`, which is session-scoped state: whether a
# given call warns depends on whether some earlier test already spent the
# budget, and therefore on file order. Force every one of them to signal so the
# suite tests the condition rather than the ordering.
options(rlib_warning_verbosity = "verbose")
