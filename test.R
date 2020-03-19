data.one.SIR <- reform.data(cevadata,
                            sampleday.vars = paste0(c(1:14),".dpch"),
                            model = "SIR",
                            positive.if = "max",
                            inf.rule = 2,
                            rec.rule = 2,
                            cut.off = 36)
data.one.SIR[[1]]
