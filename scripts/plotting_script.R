

png(file="~/Documents/Software/cicero_website/images/Cicero_example.png", width = 3, height = 2, units = "in", res=600)
plot_connections(conns, "chr18", 8575097, 8839855, gene_model = gene_annotation_sample, coaccess_cutoff = .25, connection_width = .5, collapseTranscripts = "longest" )
dev.off()


png(file="~/Documents/Software/cicero_website/images/comparison_track.png", width = 3, height = 2, units = "in", res=600)
plot_connections(conns, "chr18", 10000, 112367,
                 gene_model = gene_annotation_sample,
                 coaccess_cutoff = 0,
                 connection_width = .5,
                 comparison_track = chia_conns,
                 include_axis_track = F,
                 collapseTranscripts = "longest")
dev.off()

png(file="~/Documents/Software/cicero_website/images/viewpoint_example.png", width = 3, height = 2, units = "in", res=600)
plot_connections(conns, "chr18", 10000, 112367, viewpoint = "chr18_48000_53000",
                 gene_model = gene_annotation_sample,
                 coaccess_cutoff = 0,
                 connection_width = .5,
                 comparison_track = chia_conns,
                 include_axis_track = F,
                 collapseTranscripts = "longest")
dev.off()

png(file="~/Documents/Software/cicero_website/images/alpha_by_co_example.png", width = 3, height = 2, units = "in", res=600)
plot_connections(conns, alpha_by_coaccess = TRUE, "chr18", 8575097, 8839855, gene_model = gene_annotation_sample, coaccess_cutoff = 0.1, connection_width = .5, collapseTranscripts = "longest" )
dev.off()

png(file="~/Documents/Software/cicero_website/images/alpha_by_co_F_example.png", width = 3, height = 2, units = "in", res=600)
plot_connections(conns, alpha_by_coaccess = FALSE, "chr18", 8575097, 8839855, gene_model = gene_annotation_sample, coaccess_cutoff = 0.1, connection_width = .5, collapseTranscripts = "longest" )
dev.off()


png(file="~/Documents/Software/cicero_website/images/colors1_example.png", width = 3, height = 2, units = "in", res=600)
plot_connections(conns,
                 "chr18", 10000, 112367,
                 connection_color = "in_chia_100",
                 comparison_track = chia_conns,
                 peak_color = "green",
                 comparison_peak_color = "orange",
                 comparison_connection_color = "purple",
                 gene_model_color = "#2DD881",
                 gene_model = gene_annotation_sample,
                 coaccess_cutoff = 0.1,
                 connection_width = .5,
                 collapseTranscripts = "longest" )
dev.off()

conns$conn_color <- "orange"
conns$conn_color[conns$in_chia_100] <- "green"
png(file="~/Documents/Software/cicero_website/images/colors2_example.png", width = 3, height = 2, units = "in", res=600)
plot_connections(conns,
                 "chr18", 10000, 112367,
                 connection_color = "conn_color",
                 comparison_track = chia_conns,
                 peak_color = "green",
                 comparison_peak_color = "orange",
                 comparison_connection_color = "purple",
                 gene_model_color = "#2DD881",
                 gene_model = gene_annotation_sample,
                 coaccess_cutoff = 0.1,
                 connection_width = .5,
                 collapseTranscripts = "longest" )
dev.off()

conns$peak1_color <- "purple"
conns$peak1_color[conns$Peak1 == "chr18_11604_13986"] <- "cyan"
png(file="~/Documents/Software/cicero_website/images/colors3_example.png", width = 3, height = 2, units = "in", res=600)
plot_connections(conns,
                 "chr18", 10000, 112367,
                 connection_color = "green",
                 comparison_track = chia_conns,
                 peak_color = "peak1_color",
                 comparison_peak_color = "orange",
                 comparison_connection_color = "purple",
                 gene_model_color = "#2DD881",
                 gene_model = gene_annotation_sample,
                 coaccess_cutoff = 0.1,
                 connection_width = .5,
                 collapseTranscripts = "longest" )
dev.off()


myog_plot <- plot_connections(subb_ci, paste("chr", chr, sep=""),
                              gene_model_color = "#C7CBDB", minbp, connection_ymax = .6, comparison_ymax = .6, connection_color = "#A54657",
                              viewpoint = "chr1_203055164_203055664", viewpoint_color = "#F9DC5C", viewpoint_fill = "#F9DC5C",
                              comparison_track = subb_ci2, comparison_connection_color = "#05668D", peak_color = "#000000", comparison_peak_color = "#000000",
                              maxbp, collapseTranscripts = TRUE, include_axis_track = TRUE,
                              coaccess_cutoff = 0.25,
                              comparison_coaccess_cutoff = 0.25,
                              gene_model = gff_for_plots_sub, connection_width = .5, return_as_list = TRUE)
myog_plot@trackList <- myog_plot@trackList[c(1,4,5,3,6)]
size_track <- c(.3, .3, .25, .05, .3)

conservation <- UcscTrack(genome = "hg19", chromosome = "chr18",
                          track = "Conservation", table = "phyloP100wayAll",fontsize.group=6,fontsize=6, cex.axis=.8,
                          from = 10000, to = 112367, trackType = "DataTrack",
                          start = "start", end = "end", data = "score", size = .1,
                          type = "histogram", window = "auto", col.histogram = "#587B7F",
                          fill.histogram = "#587B7F", ylim = c(-1, 2.5),
                          name = "Conservation")


plot_list <- plot_connections(conns,
                              "chr18", 10000, 112367,
                              gene_model = gene_annotation_sample,
                              coaccess_cutoff = 0.1,
                              connection_width = .5,
                              collapseTranscripts = "longest",
                              return_as_list = TRUE)

png(file="~/Documents/Software/cicero_website/images/return_as_list_example.png", width = 3, height = 2, units = "in", res=600)
Gviz::plotTracks(plot_list,
                sizes = c(2,.5,.5,.5),
                    from = 10000, to = 112367, chromosome = "chr18",
                    transcriptAnnotation = "symbol",
                    col.axis = "black",
                    fontsize.group = 6,
                    fontcolor.legend = "black",
                    lwd=.3,
                    title.width = .5,
                    background.title = "transparent",
                    col.border.title = "transparent")

dev.off()








