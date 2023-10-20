#Read metadata
metadata <- read_excel(here::here('data-raw/tinres_limma_design.xlsx')) %>%
    dplyr::filter(!sample %in% c("s_1", "s_2", "s_3", "s_4", "s_5", "s_6", "s_7", "s_8", "s_13", "s_15", "s_50", "s_52")) %>%
    mutate_all(factor) %>%
    tibble::column_to_rownames("sample")

