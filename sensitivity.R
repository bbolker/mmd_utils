library(gamm4)
library(broom.mixed)
library(dplyr)
library(colorspace)
library(ggplot2); theme_set(theme_bw() + theme(panel.spacing = grid::unit(0, "lines")))
## library(plotly)
source("gamm4_utils.R")
ecoreg <- readRDS("ecoreg.rds")

b_ord <- (ecoreg
    %>% group_by(biome)
    %>% summarise(across(NPP_mean, mean, na.rm = TRUE), .groups="drop")
    %>% arrange(NPP_mean)
    %>% pull(biome)
    %>% as.character()
)
b_ord <- c(b_ord, "full_model")


if (!file.exists("sensitivity_fits.rds")) {
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

bestres <- list()
for (curr_taxon in names(af)) {
    ## curr_taxon <- "mamph_log"
    best_ind <- ss_best |> filter(taxon == curr_taxon) |> pull(n)
    model <- af[[curr_taxon]][[best_ind]]
    bestres[[curr_taxon]] <- (tidy(model, conf.int = TRUE, effects = "fixed")
        |> dplyr::select(term, estimate, lwr = conf.low, upr = conf.high)
        |> dplyr::filter(term != "(Intercept)")
    )
}
bestres2 <- (bind_rows(bestres, .id = "taxon")
    |> mutate(biome_ex = "full_model")
)
res2 <- bind_rows(res, .id = "taxon")

## order terms by mean estiamte
f_ord <- (res2 |> group_by(term)
    |> summarise(across(estimate, mean))
    |> arrange(estimate)
    |> ungroup()
    |> pull(term)
)


res3 <- (res2
    |> bind_rows(bestres2)
    |> mutate(across(term, factor, levels = f_ord))
    |> mutate(across(biome_ex, factor, levels = b_ord))
)



saveRDS(res3, file = "sensitivity_fits.rds")
} else {
    res3 <- readRDS("sensitivity_fits.rds")
}

green_to_tan <- c(sequential_hcl(length(b_ord)-1, h1 = 150, h2=60, c1=48, c2=63, l1=52, l2=70, rev=TRUE),
                  "#000000")
qual <- c(qualitative_hcl(n=length(b_ord)-1, rev = TRUE), "#000000")
pd <- position_dodge(width = 0.7)
gg1 <- ggplot(res3, aes(estimate, term, colour = biome_ex)) +
    geom_linerange(aes(xmin = lwr, xmax = upr), position = pd) +
    geom_point(position = pd, size = 0.7) +
    geom_vline(xintercept = 0, lty = 2) +
    ## scale_color_manual(values = green_to_tan, guide = guide_legend(reverse=TRUE)) +
    ## scale_color_discrete_qualitative(guide = guide_legend(reverse=TRUE)) +
    scale_color_manual(values = qual, guide = guide_legend(reverse=TRUE)) +
    facet_wrap(~taxon) +
    labs(y="")

print(gg1)
ggsave(gg1, file = "sensitivity.pdf", width = 12, height = 5)

if (interactive()) ggplotly()
