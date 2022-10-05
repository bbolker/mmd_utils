library(gamm4)
library(broom.mixed)
library(dplyr)
library(colorspace)
library(ggplot2); theme_set(theme_bw() + theme(panel.spacing = grid::unit(0, "lines")))
library(plotly)
source("gamm4_utils.R")
ecoreg <- readRDS("ecoreg.rds")

system.time(af <- readRDS("allfits_gamm4.rds")) ## 6-7 seconds
af <- af[1:3] ## drop plant fits

t_names <- names(af)
res <- list()

ss <- readRDS("allfits_sum_gamm4.rds")

ss_best <- (ss$sum
    |> filter(taxon != "plants_log")
    |> group_by(taxon)
    |> mutate(n = 1:n())
    |> filter(best)
    |> dplyr::select(taxon, model, n)
)

for (curr_taxon in names(af)) {
    ## curr_taxon <- "mamph_log"
    best_ind <- ss_best |> filter(taxon == curr_taxon) |> pull(n)
    model <- af[[curr_taxon]][[best_ind]]
    cc <- attr(model, "call")
    cc[[1]] <- quote(g4fit)
    ## re-run with each biome in turn eliminated, extract coefficients
    bres <- list()
    for (biome_ex in unique(ecoreg$biome)) {
        cat(curr_taxon, biome_ex, "\n")
        ## biome_ex <- "Trop/Subtrop Moist"
        rdata <- subset(ecoreg, biome != biome_ex)
        cc$data <- quote(rdata)
        new_fit <- eval(cc)
        bres[[biome_ex]] <- (tidy(new_fit, conf.int = TRUE, effects = "fixed")
            |> dplyr::select(term, estimate, lwr = conf.low, upr = conf.high)
            |> dplyr::filter(term != "(Intercept)")
        )
    }
    res[[curr_taxon]] <- dplyr::bind_rows(bres, .id = "biome_ex")
}

res2 <- (bind_rows(res, .id = "taxon")
    |> mutate(across(term, forcats::fct_reorder, estimate))
)

pd <- position_dodge(width = 0.7)
gg1 <- ggplot(res2, aes(estimate, term, colour = biome_ex)) +
    geom_linerange(aes(xmin = lwr, xmax = upr), position = pd) +
    geom_point(position = pd, size = 0.7) +
    geom_vline(xintercept = 0, lty = 2) +
    scale_color_discrete_qualitative(guide = guide_legend(reverse=TRUE)) +
    facet_wrap(~taxon)

print(gg1)
ggsave(gg1, file = "sensitivity.pdf", width = 12, height = 5)

if (interactive()) ggplotly()
