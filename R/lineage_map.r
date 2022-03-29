#' Lineage Figures - Lineage Map
#'
#' This group of Lineage Figure functions produce figures to assist in interpretation of lineage data.
#'
#' This function produces a world map showing which lineages have been found in each country.
#' The outputs of the seq_designation, node_info and lineage_info functions are required, along with the
#' phylogenetic tree and corresponding metadata file used as input for the sequence data, node info and lineage
#' information.
#'
#' @param lineage_info The output of the lineage_info function
#' @param node_data The output of the node_info function
#' @param tree A phylogenetic tree
#' @param metadata The metadata corresponding to the sequences in the tree, including "ID" "assignment" "country" and "year"
#' @param sequence_data The output of the seq_designation function
#' @param map The base layer of the map
#' @return A world map showing which lineages have been found in each country
#' @export
lineage_map <- function(lineage_info, node_data, tree, metadata, sequence_data, map) {
  tree$tip.label <- gsub("\\..*", "", tree$tip.label, perl = T)
  tree$node.comment<- gsub(".*=", "", tree$node.label, perl = T)

  lineage_info$colour<-NA

  Colours<-c("Reds","Purples","YlOrBr","PuBuGn","YlOrRd","OrRd","PuBu","Pastel1","Greens","Greys",
             "GnBu","BuGn","RdPu","Oranges","BuPu","YlGn","PuRd","YlGnBu")

  lineages<-data.frame(lineage = lineage_info$lineage, subclade = NA)

  for (i in 1:length(lineages$lineage)) {
    lineages$subclade[i]<-strsplit(lineages$lineage[i], "_")[[1]][1]
  }

  letters <- c("A1", "B1", "C1", "D1", "E1", "F1", "G1", "H1", "I1", "J1", "K1", "L1", "M1", "N1",
               "O1", "P1", "Q1", "R1", "S1", "T1", "U1", "V1", "W1", "X1", "Y1", "Z1")

  if(length(grep("_", lineage_info$lineage)) != 0) {
    if (length(which(lineages$subclade %in% letters)) != 0) {
      lineages<-lineages[-c(which(lineages$subclade %in% letters)),]
    }
  }

  clades<-unique(lineages$subclade)

  if(length(grep("\\.", clades)) != 0 ) {
    clades<-clades[-c(grep("\\.", clades))]
  }

  if (length(which(clades == "Cosmopolitan")) != 0) {
    clades<-clades[-c(which(clades == "Cosmopolitan"))]
  }

  if (length(which(clades == "Asian")) != 0) {
    clades<-clades[-c(which(clades == "Asian"))]
  }

  if (length(which(clades == "Bats")) != 0) {
    clades<-clades[-c(which(clades == "Bats"))]
  }
  lineage<-lineage_info$lineage[-c(grep("_", lineage_info$lineage))]
  cols<-RColorBrewer::brewer.pal(9, "Blues")
  pal<-colorRampPalette(c(cols))
  pal<-rev(pal(length(lineage)))
  lineage_info$colour[-c(grep("_", lineage_info$lineage))]<-pal

  for (i in 1:length(clades)) {
    lineage<-grep(clades[i], lineage_info$lineage)
    cols<-RColorBrewer::brewer.pal(3, Colours[i])
    pal<-colorRampPalette(c(cols))
    pal<-rev(pal(length(lineage)))
    lineage_info$colour[(grep(clades[i], lineage_info$lineage))]<-pal
  }

  if (map == "default"){
    world<-rgdal::readOGR("../inst/extdata/Shapefile", "world-administrative-boundaries")
  } else {
    world<-map
  }

  #' **Cleaning the data**
  #' Find which country names do not match between data and map file
  map_countries <- as.character(world@data$name)

  areas <- unique(metadata$country)
  no_match <- setdiff(areas, map_countries); message(length(no_match), " countries are mis-matched: \n", paste0(no_match, collapse="\n"))


  metadata$country[grepl("USA|United States", metadata$country)] <- "United States of America"
  metadata$country[grepl("Tanzania", metadata$country)] <- "United Republic of Tanzania"
  metadata$country[grepl("Serbia", metadata$country)] <- "Republic of Serbia"
  metadata$country[grepl("Czechia", metadata$country)] <- "Czech Republic"
  metadata$country[grepl("Ivory", metadata$country)] <- "Côte d'Ivoire" # The weird characters mean we have to do a string search!
  metadata$country[grepl("Ivoire", metadata$country)] <- "Côte d'Ivoire" # The weird characters mean we have to do a string search!
  metadata$country[grepl("Lao", metadata$country)] <- "Lao People's Democratic Republic"
  metadata$country[grepl("Laos", metadata$country)] <- "Lao People's Democratic Republic"

  areas <- unique(metadata$country)
  no_match <- setdiff(areas, map_countries); message(length(no_match), " countries are mis-matched: \n", paste0(no_match, collapse="\n"))

  no_match <- setdiff(no_match, map_countries); message(length(no_match), " countries are mis-matched: \n", paste0(no_match, collapse="\n"))

  countries <- unique(metadata$country)

  sequence_data<-sequence_data[-c(which(is.na(sequence_data$lineage))),]

  lineage_table <- data.frame(matrix(ncol = (length(unique(sequence_data$lineage))+1), nrow = length(countries)))
  x <- c("country", unique(sequence_data$lineage))
  colnames(lineage_table) <- x

  lineage_table$country<-countries

  for (i in 1:length(countries)) {
    country_place<-which(lineage_table$country == countries[i])
    lineages<-unique(sequence_data$lineage[which(sequence_data$ID %in% metadata$ID[which(metadata$country == countries[i])])])

    for (j in 1:length(lineages)) {
      country_lineage<-which(colnames(lineage_table) == lineages[j])
      lineage_table[country_place, country_lineage]<-length(which((sequence_data$lineage[which(sequence_data$ID %in%
                                                                                                 metadata$ID[which(metadata$country == countries[i])])]) == lineages[j]))
    }
  }

  lineage_table[is.na(lineage_table)] <- 0

  lineage_table$LAT<-NA
  lineage_table$LON<-NA

  if (length((which(lineage_table$country == "-"))) != 0) {
    lineage_table<-lineage_table[-c(which(lineage_table$country == "-")),]
  }

  if (length((which(lineage_table$country == "Grenada"))) != 0) {
    lineage_table<-lineage_table[-c(which(lineage_table$country == "Grenada")),]
  }

  if (length((which(lineage_table$country == "French Guiana"))) != 0) {
    lineage_table<-lineage_table[-c(which(lineage_table$country == "French Guiana")),]
  }

  for (i in 1:length(lineage_table$country)) {
    lineage_table$LAT[i]<-world@polygons[which(world@data$name == lineage_table$country[i])][[1]]@labpt[2]
    lineage_table$LON[i]<-world@polygons[which(world@data$name== lineage_table$country[i])][[1]]@labpt[1]
  }

  plot<-ggplot2::ggplot()+
    ggplot2::geom_polygon(data = world, ggplot2::aes( x = long, y = lat, group = group), fill="grey80", color="white") +
    ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())

  lineage_table <- with(lineage_table, lineage_table[abs(LON) < 150 & abs(LAT) < 70,])
  n <- nrow(lineage_table)
  lineage_table$region <- factor(1:n)
  lineage_table$radius<-NA
  for (i in 1:length(lineage_table$country)) {
    lineage_table$radius[i]<-sum(lineage_table[i,2:(ncol(lineage_table)-4)])
  }

  lineage_table$radius[which(lineage_table$radius %in% 6:10)]<-6
  lineage_table$radius[which(lineage_table$radius %in% 11:20)]<-7
  lineage_table$radius[which(lineage_table$radius %in% 21:50)]<-8
  lineage_table$radius[which(lineage_table$radius %in% 51:100)]<-9
  lineage_table$radius[which(lineage_table$radius %in% 101:200)]<-10
  lineage_table$radius[which(lineage_table$radius > 200)]<-11

  lineage_table$radius<-lineage_table$radius/1.5

  colour_table<-data.frame(lineage = colnames(lineage_table[2:(ncol(lineage_table)-4)]), colour = NA)

  for (i in 1:length(colour_table$lineage)) {
    colour_table$colour[i]<-lineage_info$colour[which(lineage_info$lineage == colour_table$lineage[i])]
  }

  plot_world<-plot + scatterpie::geom_scatterpie(ggplot2::aes(x=LON, y=LAT, group=region, r=radius),
                                     data=lineage_table, cols=c(colnames(lineage_table)[2:(ncol(lineage_table)-4)]), color=NA, alpha=.8)+
    ggplot2::theme(legend.position = "none") + ggplot2::scale_fill_manual(values = c(colour_table$colour))+
    ggplot2::theme(panel.grid.major = ggplot2::element_blank())+
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          rect = ggplot2::element_blank(),
          axis.title.y=ggplot2::element_blank(),
          axis.title.x=ggplot2::element_blank())


  return(plot_world)
}


