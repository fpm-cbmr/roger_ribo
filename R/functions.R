
#' Volcano
#'
#' @param dataset
#' @param type (p-value adjustment. Either p, xiao, or q)
#' @param threshold (cut-off)
#'
#' @return a volcano plot

volcano <- function (dataset, type, threshold) {

    dataset%>%
        mutate(color = case_when(
            logFC >= 0 & {{type}} <= {{threshold}} ~ "Upregulated",
            logFC <= 0 & {{type}} <= {{threshold}} ~ "Downregulated",
            TRUE ~ "Unchanged")) %>%
        ggplot(aes(x=logFC, y=-log10(P.Value), label=row.names(dataset)))+
        geom_point(aes(color = color, alpha=color), size = 3)+
        theme(panel.background = element_rect(fill = "transparent", colour = NA),
              plot.background = element_rect(fill = "transparent", colour = NA),
              panel.grid.minor=element_blank(),
              panel.grid.major = element_blank(),
              axis.line = element_line(colour = "black"),
              legend.text=element_text(size=10),
              text = element_text(size = 12),
              legend.title = element_blank(),
              legend.key = element_blank(),
              plot.title = element_text(hjust = 0.5, face="bold"),
              axis.text.x = element_text(color="black", size=10),
              axis.text.y = element_text(color="black", size=10),
              legend.position = "none")+
        geom_text_repel(point.size=4, size=3, min.segment.length = 0.1, force=0.3)+
        scale_color_manual(breaks = c("Upregulated", "Downregulated", "Unchanged"),
                           values=c("dodgerblue3", "firebrick3", "gray50"))+
        scale_alpha_manual(breaks = c("Upregulated", "Downregulated", "Unchanged"),
                           values=c(1, 1, 0.1))+
        xlab("Log2fold difference (type IIa - type I)") + ylab("-log10(p)")
}


#' Volcano for interaction plots
#'
#' @param dataset
#' @param type (p-value adjustment. Either p, xiao, or q)
#' @param threshold (cut-off)
#'
#' @return a volcano plot

volcano_interaction <- function (dataset, type, threshold) {

    dataset%>%
        mutate(color = case_when(
            logFC >= 0 & {{type}} <= {{threshold}} ~ "Upregulated",
            logFC <= 0 & {{type}} <= {{threshold}} ~ "Downregulated",
            TRUE ~ "Unchanged")) %>%
        ggplot(aes(x=logFC, y=-log10(P.Value), label=row.names(dataset)))+
        geom_point(aes(color = color, alpha=color), size = 3)+
        theme(panel.background = element_rect(fill = "transparent", colour = NA),
              plot.background = element_rect(fill = "transparent", colour = NA),
              panel.grid.minor=element_blank(),
              panel.grid.major = element_blank(),
              axis.line = element_line(colour = "black"),
              legend.text=element_text(size=10),
              text = element_text(size = 12),
              legend.title = element_blank(),
              legend.key = element_blank(),
              plot.title = element_text(hjust = 0.5, face="bold"),
              axis.text.x = element_text(color="black", size=10),
              axis.text.y = element_text(color="black", size=10),
              legend.position = "none")+
        geom_text_repel(point.size=4, size=3, min.segment.length = 0.1, force=0.3)+
        scale_color_manual(breaks = c("Upregulated", "Downregulated", "Unchanged"),
                           values=c("dodgerblue3", "firebrick3", "gray50"))+
        scale_alpha_manual(breaks = c("Upregulated", "Downregulated", "Unchanged"),
                           values=c(1, 1, 0.1))+
        xlab("Log2 difference (type II - type I)") + ylab("-log10(p)")+
        xlim(-2.5, 2.5)+
        ylim(0,9)
}



#' Volcano with continuous color based on type
#'
#' @param dataset
#' @param type
#' @param threshold
#'
#' @return a volcano plot

volcano_con <- function (dataset, type) {

    xiao_gradient <- colorRampPalette(c((colorRampPalette(c("gray","gray", "gray", "gray","gray", "#F94040"))(50)), rev(colorRampPalette(c("gray","gray","gray","gray","gray", "#5757F9"))(50))))

    dataset%>%
        dplyr::mutate(color = ifelse(logFC>0, {{type}}, {{type}}*-1)) %>%
        ggplot(aes(x=logFC, y=-log10(P.Value), label=row.names(dataset)))+
        geom_point(aes(color = color), size = 3)+
        theme(panel.background = element_rect(fill = "transparent", colour = NA),
              plot.background = element_rect(fill = "transparent", colour = NA),
              panel.grid.minor=element_blank(),
              panel.grid.major = element_blank(),
              axis.line = element_line(colour = "black"),
              legend.text=element_text(size=10),
              text = element_text(size = 12),
              legend.title = element_blank(),
              legend.key = element_blank(),
              plot.title = element_text(hjust = 0.5, face="bold"),
              axis.text.x = element_text(color="black", size=10),
              axis.text.y = element_text(color="black", size=10))+
        geom_text_repel(point.size=4, size=3, min.segment.length = 0.1, force=0.3)+
        xlab("Log2fold change") + ylab("-log10(p)")+
        scale_colour_gradientn(colors=xiao_gradient(100))
}


#' Title
#'
#' @param se
#'
#' @return Returns a heatmap visualizing missing values (in white) and valid values (in black) and column annotations from metadata in the form of intervention, time, id, and fiber type.
#' @export Nothing
#'
#' @examples missing_plot(se_ter_i)
missing_plot <- function(se){

#Create a new data frame containing only proteins with missing values
df_missing <- SummarizedExperiment::assay({{se}}) %>%
    dplyr::filter(!complete.cases(.)) %>% #Removes rows with no missing values
    dplyr::mutate_all(~ ifelse(is.na(.), 0, 1)) #NA's replaced with 0 and everything else with 1

#Heatmap
heatmap_missing <- pheatmap(df_missing,
                            cluster_rows = F,
                            cluster_cols = T,
                            annotation_col = dplyr::select(metadata, c("fiber_type", "time", "intervention", "id")),
                            annotation_colors=list(fiber_type=c(I="#db464b",II="#7798cf"),
                                                   time=c(pre="#ff6361", post="#ffa600"),
                                                   intervention=c(terbutaline="#c06c85", resistance="#345c7e")),
                            show_rownames = F,
                            color=colorRampPalette(c("white", "black"))(2),
                            #cellwidth =5,
                            border_color = NA,
                            legend = F,
                            annotation_legend = T
)
}


#' Title
#'
#' @param summarized_df e.g. df_long_l2fc_mean (string)
#' @param intervention e.g. terbutaline or resistance (string)
#' @param fiber_type e.g. I or II (string)
#'
#' @return a bar plot of mitochondrial subunits ordered by median l2fc.
#' @export
#'
#' @examples
mito_subunits_median <- function(summarized_df, intervention, fiber_type) {

df_long_l2fc_mean %>%
    dplyr::filter(intervention == {{intervention}} & fiber_type == {{fiber_type}}) %>%
    dplyr::filter(grepl("subunit", mito, ignore.case = T)) %>%
    dplyr::filter(!is.na(l2fc_median)) %>%
    ggplot2::ggplot(aes(x=reorder(protein, l2fc_median), y=l2fc_median, fill=fiber_type))+
    geom_col()+
    theme(
        axis.text.x= element_text(angle = 90, vjust = 0.5, hjust=1, color="black", size=10),
        axis.text.y= element_text(color="black", size=10),
        axis.line = element_line(color="black", size = 0.5),
        axis.ticks = element_blank(),
        panel.background = element_rect(color = "black", fill=NA, size = 0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        plot.title = element_text(size = 12),
        text = element_text(family = "Source Sans Pro", size=11),
        legend.position = "none"
    )+
    scale_fill_manual(values = c("I" = "#cb4c52",
                                 "II" = "#7c98ce"))+
    labs(y = "Median l2fc", x = "")+
    ylim(-.8, 0.3)
}



#' Title
#'
#' @param summarized_df
#'
#' @return A violion plot of median l2fc of mitochondrial subunits stacked.
#' @export
#'
#' @examples
mito_subunits_summed_median <- function (summarized_df) {

    df_long_l2fc_mean %>%
    dplyr::filter(grepl("subunit", mito, ignore.case = T)) %>%
    ggplot2::ggplot(aes(x=intervention, y=l2fc_median, fill=fiber_type))+
    geom_violin(aes(group=interaction(fiber_type, intervention)), position =position_dodge(width=1))+
    geom_point(aes(fill=fiber_type), na.rm=TRUE, position = position_dodge(width = 1), size=5, alpha=0.2)+
    geom_hline(yintercept = 0, linetype = "dashed")+
    theme(
        axis.text.x= element_text(color="black", size = 12),
        axis.text.y= element_text(color="black", size = 12),
        axis.line = element_line(color="black", size = 0.5),
        panel.background = element_rect(color = "black", fill=NA, size = 0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        plot.title = element_text(size = 12),
        plot.margin = unit(c(5, 5, 5, 25), "pt"),
        text = element_text(family = "Source Sans Pro", size=12)
    )+
    scale_fill_manual(values=c("#cb4c52", "#7c98ce"),
                      name = "Fiber Type",
                      labels = c("I",
                                 "II"))+
    ggtitle("Mitochondrial subunits")+
    labs(x = "", y = "Summed median l2fc")+
    scale_x_discrete(labels = c(resistance = "Resistance training", terbutaline = expression("Beta"[2]*"-agonist")))
}

#' Title
#'
#' @param summarized_df e.g. df_long_l2fc_mean (string)
#' @param intervention e.g. terbutaline or resistance (string)
#' @param fiber_type e.g. I or II (string)
#'
#' @return a bar plot of mitochondrial subunits ordered by median l2fc.
#' @export
#'
#' @examples
ribo_subunits_median <- function(summarized_df, intervention, fiber_type) {

    df_long_l2fc_mean %>%
        dplyr::filter(intervention == {{intervention}} & fiber_type == {{fiber_type}}) %>%
        dplyr::filter(grepl('cytosolic ribosome', gocc, ignore.case = T)) %>%
        dplyr::filter(!is.na(l2fc_median)) %>%
        ggplot2::ggplot(aes(x=reorder(protein, l2fc_median), y=l2fc_median, fill=fiber_type))+
        geom_col()+
        theme(
            axis.text.x= element_text(angle = 90, vjust = 0.5, hjust=1, color="black", size=10),
            axis.text.y= element_text(color="black", size=10),
            axis.line = element_line(color="black", size = 0.5),
            axis.ticks = element_blank(),
            panel.background = element_rect(color = "black", fill=NA, size = 0.5),
            panel.grid.minor=element_blank(),
            panel.grid.major = element_blank(),
            plot.background = element_blank(),
            plot.title = element_text(size = 12),
            text = element_text(family = "Source Sans Pro", size=11),
            legend.position = "none"
        )+
        scale_fill_manual(values = c("I" = "#cb4c52",
                                     "II" = "#7c98ce"))+
        labs(y = "Median l2fc", x = "")+
        ylim(-.5, 0.65)
}

#' Title
#'
#' @param summarized_df
#'
#' @return A violion plot of median l2fc of mitochondrial subunits stacked.
#' @export
#'
#' @examples
ribo_subunits_summed_median <- function (summarized_df) {

    df_long_l2fc_mean %>%
        dplyr::filter(grepl('cytosolic ribosome', gocc, ignore.case = T)) %>%
        ggplot2::ggplot(aes(x=intervention, y=l2fc_median, fill=fiber_type))+
        geom_violin(aes(group=interaction(fiber_type, intervention)), position =position_dodge(width=1))+
        geom_point(aes(fill=fiber_type), na.rm=TRUE, position = position_dodge(width = 1), size=5, alpha=0.2)+
        geom_hline(yintercept = 0, linetype = "dashed")+
        theme(
            axis.text.x= element_text(color="black", size = 12),
            axis.text.y= element_text(color="black", size = 12),
            axis.line = element_line(color="black", size = 0.5),
            panel.background = element_rect(color = "black", fill=NA, size = 0.5),
            panel.grid.minor=element_blank(),
            panel.grid.major = element_blank(),
            plot.background = element_blank(),
            plot.title = element_text(size = 12),
            plot.margin = unit(c(5, 5, 5, 25), "pt"),
            text = element_text(family = "Source Sans Pro", size=12)
        )+
        scale_fill_manual(values=c("#cb4c52", "#7c98ce"),
                          name = "Fiber Type",
                          labels = c("I",
                                     "II"))+
        ggtitle("Ribosomal proteins")+
        labs(x = "", y = "Summed median l2fc")+
        scale_x_discrete(labels = c(resistance = "Resistance training", terbutaline = expression("Beta"[2]*"-agonist")))
}
