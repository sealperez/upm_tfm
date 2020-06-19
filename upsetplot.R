## Remember to set your working directory!!
## setwd("path_to_your_WD")

library('dplyr')
library('ComplexHeatmap')

gene_presence_absence <- read.csv("gene_presence_absence.csv", encoding="437", row.names=1, na.strings="")
df<-gene_presence_absence[,14:19]
for(i in 1:ncol(df)){
	if(is.factor(df[,i]))
		levels(df[,i]) <- c(levels(df[,i]),"0","1")
}
df[!is.na(df)]<-"1"
df[is.na(df)] <- "0"
genes<-rownames(df)
df<-apply(df[,1:ncol(df)],2, function(x) as.numeric(as.character(x)))
df<-as.data.frame(df)
rownames(df)<-genes
names(df)<-c('R_australiborealis_S264','R_collisarenosi_S260','R_epipactidis_S256','R_gaditana_S61','R_nectarea_B1A','R_nectarea_8N4')

m<-make_comb_mat(df)

pdf("upsetplot.pdf", width = 16, height = 9)

ss = set_size(m)
cs = comb_size(m)
ht = UpSet(m, 
    set_order = order(ss),
    comb_order = order(-cs),
    bg_col = c("#F0F0FF", "#FFF0F0"), bg_pt_col = "#CCCCFF",
    top_annotation = HeatmapAnnotation(
        "Gene intersections" = anno_barplot(cs, 
            ylim = c(0, max(cs)*1.1),
            border = FALSE, 
            gp = gpar(fill = "steelblue"), 
            height = unit(4, "cm")
        ), 
        annotation_name_side = "left", 
        annotation_name_rot = 90),
    left_annotation = rowAnnotation(
        "No. genes" = anno_barplot(-ss, 
            baseline = 0,
            axis_param = list(
                at = c(0, -1000, -2000, -3000),
                labels = c(0, 1000, 2000, 3000),
                labels_rot = 0),
            border = FALSE, 
            gp = gpar(fill = "lightgreen"), 
            width = unit(4, "cm")
        ),
        set_name = anno_text(set_name(m), 
            location = 0.5, 
            just = "center",
            width = max_text_width(set_name(m)) + unit(4, "mm"))
    ), 
    right_annotation = NULL,
    show_row_names = FALSE)
ht = draw(ht)
od = column_order(ht)
decorate_annotation("Gene intersections", {
    grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
        default.units = "native", just = c("left", "bottom"), 
        gp = gpar(fontsize = 10, col = "black"), rot = 45)
})

dev.off()
