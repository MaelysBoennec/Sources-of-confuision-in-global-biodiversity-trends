################################################################################
### Maelys Boennec
### Revised in 04/2023
################################################################################

###############################################################################
###=========================== LOAD LIBRARIES ==============================###
###############################################################################

library(ggplot2)
library(tidyr)
library("gridExtra")
library(cowplot)
library(ggalluvial)
library(forcats)
library(ggfittext)
library(ggbreak)
library(likert)
library(ggplotify)
library(ggfortify)
library(RColorBrewer)
library(devtools)
library(moonBook)
library(webr)


###############################################################################
###============================== LOAD DATA ================================###
###############################################################################


# Load meta data ----------------------------------------------------------

stats <- readr::read_csv2(
  here::here("data","metadata_dup_v2.csv")
)

# rename some columns
stats$temporal_range <- stats$temporal_extent
stats$taxonomic_coverage <- stats$taxonomic_groups
stats$metric_method <- stats$method

# charge an older version
# only some categories name differ but this was used for figure 2 initially
# we didn't change the script to produce the figure with the actualized metadata version
metadata <- readr::read_csv2(
  here::here("data","metadata_dup.csv")
)

# Categorisation process --------------------------------------------------

## taxonomic categories
stats$n_taxo <- sapply(strsplit(stats$taxonomic_coverage, ", "), length)
stats$taxo <- ifelse(stats$n_taxo==1,"1 group",NA)
stats$taxo <- ifelse(stats$n_taxo>1 & stats$n_taxo<4,"2-3 groups",stats$taxo)
stats$taxo <- ifelse(stats$n_taxo>=4,"4 + groups",stats$taxo)
stats$taxo <- ifelse(stats$taxonomic_coverage=="/",NA,stats$taxo)
stats$taxo <- ifelse(stats$taxonomic_coverage=="unknown","group unknown",stats$taxo)

## temporal coverage
temp <- as.data.frame(stringr::str_split(stats$temporal_range,"-",simplify=TRUE))
stats$start_year <- as.numeric(as.character(temp$V1))
stats$end_year <- as.numeric(as.character(temp$V2))
stats$time_period <- stats$end_year - stats$start_year

## temporal categories
stats$time <- ifelse(stats$time_period<31,"- 30 years",NA)
stats$time <- ifelse(stats$time_period<51 & stats$time_period>30,"31 - 50 years",stats$time)
stats$time <- ifelse(stats$time_period>50,"51 + years",stats$time)
stats$time <- ifelse(stats$temporal_range=="unknown","time unknown",stats$time)

## clear conclusions
stats <- stats %>%
  dplyr::mutate(ccl2=ifelse(
    simple_conclusion=="factor-dependent",
    stringr::str_c(simple_conclusion,sub_type,sep=" "),
    as.character(simple_conclusion)))

metadata <- metadata %>%
  dplyr::mutate(ccl2=ifelse(
    simple_conclusion=="it depends",
    stringr::str_c(simple_conclusion,sub_type,sep=" "),
    as.character(simple_conclusion)))

# Subsets selection -------------------------------------------------------

# Biodiversity analysis
biodiv <- stats %>%
  dplyr::filter(type=="Empirical analysis") # Empirical biodiversity analysis only

# Number of species
biodiv$n_sp <- as.numeric(as.character(biodiv$n_species))
biodiv$sp <- ifelse(biodiv$n_sp<100,"- 100 species",NA)
biodiv$sp <- ifelse(biodiv$n_sp>=100 & biodiv$n_sp<1000,"100-1000 species",biodiv$sp)
biodiv$sp <- ifelse(biodiv$n_sp>=1000 & biodiv$n_sp<10000,"1000-10000 species",biodiv$sp)
biodiv$sp <- ifelse(biodiv$n_sp>=10000,"10000 + species",biodiv$sp)
biodiv$sp[which(biodiv$n_species=="unknown")] <- "unknown"

# Create a dataframe for figures 4B and 5B
methest <- as.data.frame(table(biodiv$data,
                               biodiv$sub_type,
                               biodiv$method,
                               biodiv$ecological_level,
                               biodiv$ccl2)) %>%
  dplyr::filter(Freq>0) %>%
  dplyr::rename(data="Var1",
                analysis="Var2",
                method="Var3",
                level="Var4",
                ccl="Var5",
                n="Freq") %>%
  dplyr::mutate(data=factor(data,levels=c("Biotime","LPD","Biotime, LPD","GPDD","Other","Aggregation")),
                method=factor(method,levels=c("linear model","global indicator","other","CC","anthropogenic drivers","anthropogenic","conservation")),
                level=factor(level,levels=c("population","species","grouped")),
                ccl=factor(ccl,levels=c("decreasing trends","factor-dependent trend","mixed trends","increasing trends",
                                        "negative effect of drivers","factor-dependent drivers","none/positive effect of drivers"))
  )

###############################################################################
###========================== GLOBAL OVERVIEW ==============================###
###############################################################################

# Number of articles per publication year ---------------------------------

## All articles together
count_year <- stats %>%
  dplyr::count(year) %>%
  # add years with 0 articles
  dplyr::right_join(y=data.frame(all_years=seq(1991,2022,1)),by = c("year" = "all_years")) %>%
  # replace NAs by 0s
  dplyr::mutate(n=ifelse(is.na(n),0,n)) %>%
  dplyr::arrange(year)

# cumulative number (multiplication for plotting purposes)
count_year$n_cum <- cumsum(count_year$n)*(14/97)

stats$type2<-as.factor(stats$type)
count_type <- data.frame(table(stats$type2,stats$year)) %>%
  dplyr::rename(type="Var1",
                year="Var2",
                n="Freq")

count_type$year <- as.numeric(as.character(count_type$year)) # transform col types

l_years <- seq(1991,2022,1) # list of all years between 1991 and 2022
l_type <- c("Empirical analysis","Review","Report","Methodological paper") # list of articles types

# add years with 0 articles for each article type
for (i in l_years){
  if (!(i %in% count_type$year)){
    for (t in l_type){
      count_type <- rbind(count_type,c(t,i,0))
    }
  }
}

# transform col types
count_type$year <- as.numeric(as.character(count_type$year))
count_type$n <- as.numeric(as.character(count_type$n))
levels(count_type$type) <- c("Biodiversity empirical analysis", "Methodological papers", "Reports", "Reviews")

## Plot the number of articles regarding the publication year
cols <- hcl.colors(4, "Geyser")
(plot_publication_year <- ggplot() +
  geom_area(data=count_type, aes(x=year, y=n, fill=type)) +
  geom_line(data=count_year, aes(x=year,y=n_cum,colour="Cumulative number of articles published"),linetype="dashed") +
  ylab("Number of articles published")+xlab("Publication year") +
  scale_y_continuous(breaks=seq(0,15,2),sec.axis = sec_axis(~.*(97/14), name="Cumulative number of articles published")) +
  scale_x_continuous(breaks=seq(1991,2022,2)) +
  scale_color_manual(values=c("grey50")) +
  scale_fill_manual(values=cols) +
  geom_line(data=count_year, aes(x=year,y=n)) +
  theme_classic() +
  theme(axis.text.x=element_text(size=10,angle=90,hjust=1,family="Montserrat"),
        axis.text.y=element_text(size=12,family="Montserrat"),
        axis.line.y.right = element_line(),
        legend.title=element_blank(),
        legend.text=element_text(size=12,family="Montserrat",face="italic"),
        legend.position = c(0.3,0.8),
        axis.title = element_text(size=12,family="Montserrat",face="bold")))

# Temporal coverage -------------------------------------------------------

cols <- hcl.colors(7, "Geyser")

## Temporal representation
biodiv_plot <- drop_na(biodiv,start_year,end_year)
biodiv$data <-as.factor(biodiv$data)

biodiv$sort <- biodiv$data
biodiv$sort <- factor(biodiv$sort, levels = c("Aggregation","Biotime","Biotime, LPD","LPD","GPDD","Other"),
                                  labels = c(1, 2, 3, 4, 5, 6))

biodiv$sort <- paste0(biodiv$sort, biodiv$start_year)
biodiv$sort <- as.numeric(as.character(biodiv$sort))

biodiv_plot<-biodiv %>% dplyr::filter(!is.na(start_year)&start_year>1800)
biodiv_plot$id <- as.character(1:nrow(biodiv_plot))

(timeline <- ggplot() +
    geom_linerange(data = biodiv_plot, aes(ymin = start_year, ymax = end_year, colour = data,
                                          x = fct_reorder(id, dplyr::desc(sort))),
                   size = 1) +
    scale_colour_manual(values = cols) +
    theme_classic() +
    coord_flip() +
    labs(x = NULL, y = "Year") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(),
          legend.text = element_text(size=12,family="Montserrat",face="italic"),
          panel.border = element_blank(),
          legend.title = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.text = element_text(size=12,family="Montserrat"),
          axis.title = element_text(size = 12,family="Montserrat",face="bold")))


# Groups coverage ---------------------------------------------------------

# create a table w/ the number of articles for each possible combination
t <- as.data.frame(table(biodiv$taxonomic_coverage)) %>%
  dplyr::filter(Freq>0)

# listing all groups considered
l_groups <- c("mammals","birds","amphibians","reptiles","fishes","insects","invertebrates","plants","unknown")

# initialize a data frame which will contain how many articles concern
# each group alone or combined
df_groups <- data.frame("group"=l_groups,
                  "alone"=rep(0,length(l_groups)),
                  "combined"=rep(0,length(l_groups)),
                  "perc"=rep(0,length(l_groups)))

# separate the lists into vectors
t$Var2 <- sapply(t$Var1,FUN=function(x){stringr::str_split(x,", ",simplify=TRUE)})

# articles considering only one group
df_groups$alone <- (dplyr::left_join(df_groups,t,by = c("group" = "Var1")))$Freq

# delete them
t$n_group <- sapply(t$Var2, length)
new_t <- t %>% dplyr::filter(n_group>1)

# articles considering several groups
df_groups$combined <- sapply(l_groups,FUN=function(y){sum(sapply(new_t$Var2,FUN=function(x){sum(y %in% x)})*new_t$Freq)})

# calculate the percentage of representation of each group regarding the 52 articles in total
df_groups$class <- c(rep("vertebrates",5),rep("invertebrates",3),"unknown")
df_groups$tot <- df_groups$alone+df_groups$combined

data <- df_groups %>%
  gather(key="observation", value="value", -c(1,4,5)) %>%
  dplyr::arrange(value)

data$group <- fct_relevel(data$group,rev(c("unknown","mammals","birds","amphibians",
                                       "reptiles","fishes",
                                         "plants","insects","invertebrates")))

data$observation <- fct_relevel(data$observation,c("combined","alone"))

# Load icons
reptiles <- png::readPNG(here::here("data","snake.png"))
plants <- png::readPNG(here::here("data","plant.png"))
mammals <- png::readPNG(here::here("data","tiger.png"))
invertebrates <- png::readPNG(here::here("data","jellyfish.png"))
insects <- png::readPNG(here::here("data","insect.png"))
fishes <- png::readPNG(here::here("data","fish.png"))
birds <- png::readPNG(here::here("data","bird.png"))
amphibians <- png::readPNG(here::here("data","frog.png"))
reptiles <- grid::rasterGrob(reptiles, interpolate = TRUE)
plants <- grid::rasterGrob(plants, interpolate = TRUE)
mammals <- grid::rasterGrob(mammals, interpolate = TRUE)
invertebrates <- grid::rasterGrob(invertebrates, interpolate = TRUE)
insects <- grid::rasterGrob(insects, interpolate = TRUE)
fishes <- grid::rasterGrob(fishes, interpolate = TRUE)
birds <- grid::rasterGrob(birds, interpolate = TRUE)
amphibians <- grid::rasterGrob(amphibians, interpolate = TRUE)

# Plot
cols <- hcl.colors(4, "Geyser")
(plot_groups <- ggplot() +
  geom_col(data=subset(data,group!="unknown" & observation!="tot"), aes(x = group, y = value, fill = observation)) +
  coord_flip() +
  scale_fill_manual(values=c(rev(cols)[3],rev(cols[4])))+
  geom_hline(yintercept=0) + scale_y_continuous(breaks=seq(0,30,5))+
  xlab("") + ylab("Number of articles") +
  expand_limits(y=c(0,32)) +
  geom_text(data=subset(df_groups, group!="unknown") ,aes(y=tot, x = group, label = group), hjust =1.1, colour = "black", family = "Montserrat",size = 4,fontface="bold") +
  annotation_custom(mammals, xmin=6.5, xmax=9.5, ymin=df_groups[which(df_groups$group=="mammals"),]$tot, ymax=df_groups[which(df_groups$group=="mammals"),]$tot+3) +
  annotation_custom(birds, xmin=5.5, xmax=8.5, ymin=df_groups[which(df_groups$group=="birds"),]$tot, ymax=df_groups[which(df_groups$group=="birds"),]$tot+3) +
  annotation_custom(amphibians, xmin=4.5, xmax=7.5, ymin=df_groups[which(df_groups$group=="amphibians"),]$tot, ymax=df_groups[which(df_groups$group=="amphibians"),]$tot+3) +
  annotation_custom(reptiles, xmin=3.5, xmax=6.5, ymin=df_groups[which(df_groups$group=="reptiles"),]$tot, ymax=df_groups[which(df_groups$group=="reptiles"),]$tot+3) +
  annotation_custom(fishes, xmin=2.5, xmax=5.5, ymin=df_groups[which(df_groups$group=="fishes"),]$tot, ymax=df_groups[which(df_groups$group=="fishes"),]$tot+3) +
  annotation_custom(plants, xmin=1.5, xmax=4.5, ymin=df_groups[which(df_groups$group=="plants"),]$tot, ymax=df_groups[which(df_groups$group=="plants"),]$tot+3) +
  annotation_custom(insects, xmin=0.5, xmax=3.5, ymin=df_groups[which(df_groups$group=="insects"),]$tot, ymax=df_groups[which(df_groups$group=="insects"),]$tot+3) +
  annotation_custom(invertebrates, xmin=-0.5, xmax=2.5, ymin=df_groups[which(df_groups$group=="invertebrates"),]$tot, ymax=df_groups[which(df_groups$group=="invertebrates"),]$tot+3) +

  theme_classic() +
  theme(legend.position = c(0.9,0.1),
        legend.title = element_blank(),
        legend.text = element_text(size=12,family="Montserrat",face="italic"),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x.top = element_blank(),
        axis.line.x.bottom = element_line(),
        axis.ticks.x.bottom = element_line(color="black"),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.text.x.bottom =element_text(size=12,family="Montserrat"),
        axis.text.y=element_blank(),
        axis.title.x = element_text(size=12,family="Montserrat",face="bold")))


# Ecological level representation -----------------------------------------

(ecolvl<-ggplot(data = methest,
                aes(axis1 = level, y=n)) +
    stat_alluvium(aes(fill = level), alpha=.7, lode.guidance = "forward") +
    geom_stratum(alpha=.25) +
    scale_fill_manual(values=rev(cols),na.value = "white") +

    annotate ("text",x = -1,y = 40,label = "Population",size=4,family="Montserrat",fontface="bold") +
    annotate ("text",x = -1,y = 18,label = "Species",size=4,family="Montserrat",fontface="bold") +
    annotate ("text",x = -1,y = 5,label = "Grouped",size=4,family="Montserrat",fontface="bold") +

    annotate ("text",x = -0.95,y = 40,label = "50% (n=24)",size=4,family="Montserrat",fontface="italic") +
    annotate ("text",x = -0.95,y = 18,label = "31% (n=15)",size=4,family="Montserrat",fontface="italic") +
    annotate ("text",x = -0.95,y = 5,label = "19% (n=9)",size=4,family="Montserrat",fontface="italic") +

    # visualization parameters
    coord_flip() + scale_y_reverse() + scale_x_reverse() +
    theme_void() +
    theme(legend.position = "none"))


# Figure 2 ----------------------------------------------------------------

ab <- plot_grid(timeline,plot_groups, labels=c("A","B"), ncol = 2, nrow = 1)
(fig2 <- plot_grid(ab,ecolvl,labels = c("","C"),ncol=1,nrow=2,rel_heights = c(3, 0.5)))


###############################################################################
###========================= IMPACT OF APPROACH ============================###
###############################################################################


# Define pie donut function -----------------------------------------------

PieDonutCustom <- function (data, mapping, start = getOption("PieDonut.start",
                                                             0), addPieLabel = TRUE, addDonutLabel = TRUE, showRatioDonut = TRUE,
                            showRatioPie = TRUE, ratioByGroup = TRUE, showRatioThreshold = getOption("PieDonut.showRatioThreshold",
                                                                                                     0.02), labelposition = getOption("PieDonut.labelposition",
                                                                                                                                      2), labelpositionThreshold = 0.1, r0 = getOption("PieDonut.r0",
                                                                                                                                                                                       0.3), r1 = getOption("PieDonut.r1", 1), r2 = getOption("PieDonut.r2",
                                                                                                                                                                                                                                              1.2), explode = NULL, selected = NULL, explodePos = 0.1,
                            color = "white", pieAlpha = 0.8, donutAlpha = 1, maxx = NULL,
                            showPieName = TRUE, showDonutName = FALSE, title = NULL,
                            pieLabelSize = 4, donutLabelSize = 3, titlesize = 5, explodePie = TRUE,
                            explodeDonut = FALSE, use.label = TRUE, use.labels = TRUE,
                            family = getOption("PieDonut.family", ""), palette_name="Dark2")
{
  (cols = colnames(data))
  if (use.labels)
    data = moonBook::addLabelDf(data, mapping)
  count <- NULL
  if ("count" %in% names(mapping))
    count <- moonBook::getMapping(mapping, "count")
  count
  pies <- donuts <- NULL
  (pies = moonBook::getMapping(mapping, "pies"))
  if (is.null(pies))
    (pies = moonBook::getMapping(mapping, "pie"))
  if (is.null(pies))
    (pies = moonBook::getMapping(mapping, "x"))
  (donuts = moonBook::getMapping(mapping, "donuts"))
  if (is.null(donuts))
    (donuts = moonBook::getMapping(mapping, "donut"))
  if (is.null(donuts))
    (donuts = moonBook::getMapping(mapping, "y"))
  if (!is.null(count)) {
    df <- data %>% group_by(.data[[pies]]) %>% dplyr::summarize(Freq = sum(.data[[count]]))
    df
  }
  else {
    df = data.frame(table(data[[pies]]))
  }
  colnames(df)[1] = pies
  df$end = cumsum(df$Freq)
  df$start = dplyr::lag(df$end)
  df$start[1] = 0
  total = sum(df$Freq)
  df$start1 = df$start * 2 * pi/total
  df$end1 = df$end * 2 * pi/total
  df$start1 = df$start1 + start
  df$end1 = df$end1 + start
  df$focus = 0
  if (explodePie)
    df$focus[explode] = explodePos
  df$mid = (df$start1 + df$end1)/2
  df$x = ifelse(df$focus == 0, 0, df$focus * sin(df$mid))
  df$y = ifelse(df$focus == 0, 0, df$focus * cos(df$mid))
  df$label = df[[pies]]
  df$ratio = df$Freq/sum(df$Freq)
  if (showRatioPie) {
    df$label = ifelse(df$ratio >= showRatioThreshold, paste0(df$label,
                                                             "\n(", scales::percent(df$ratio), ")"),
                      as.character(df$label))
  }
  df$labelx = (r0 + r1)/2 * sin(df$mid) + df$x
  df$labely = (r0 + r1)/2 * cos(df$mid) + df$y
  if (!is.factor(df[[pies]]))
    df[[pies]] <- factor(df[[pies]])
  df
  mainCol = RColorBrewer::brewer.pal(nrow(df), name=palette_name)
  df$radius = r1
  df$radius[df$focus != 0] = df$radius[df$focus != 0] + df$focus[df$focus !=
                                                                   0]
  df$hjust = ifelse((df$mid%%(2 * pi)) > pi, 1, 0)
  df$vjust = ifelse(((df$mid%%(2 * pi)) < (pi/2)) | (df$mid%%(2 *
                                                                pi) > (pi * 3/2)), 0, 1)
  df$segx = df$radius * sin(df$mid)
  df$segy = df$radius * cos(df$mid)
  df$segxend = (df$radius + 0.05) * sin(df$mid)
  df$segyend = (df$radius + 0.05) * cos(df$mid)
  df
  if (!is.null(donuts)) {
    subColor = makeSubColor(mainCol, no = length(unique(data[[donuts]])))
    subColor
    data
    if (!is.null(count)) {
      df3 <- as.data.frame(data[c(donuts, pies, count)])
      colnames(df3) = c("donut", "pie", "Freq")
      df3
      df3 <- eval(parse(text = "complete(df3,donut,pie)"))
      df3$Freq[is.na(df3$Freq)] = 0
      if (!is.factor(df3[[1]]))
        df3[[1]] = factor(df3[[1]])
      if (!is.factor(df3[[2]]))
        df3[[2]] = factor(df3[[2]])
      df3 <- df3 %>% arrange(.data$pie, .data$donut)
      a <- df3 %>% spread(.data$pie, value = .data$Freq)
      a = as.data.frame(a)
      a
      rownames(a) = a[[1]]
      a = a[-1]
      a
      colnames(df3)[1:2] = c(donuts, pies)
    }
    else {
      df3 = data.frame(table(data[[donuts]], data[[pies]]),
                       stringsAsFactors = FALSE)
      colnames(df3)[1:2] = c(donuts, pies)
      a = table(data[[donuts]], data[[pies]])
      a
    }
    a
    df3
    df3$group = rep(colSums(a), each = nrow(a))
    df3$pie = rep(1:ncol(a), each = nrow(a))
    total = sum(df3$Freq)
    total
    df3$ratio1 = df3$Freq/total
    df3
    if (ratioByGroup) {
      df3$ratio = scales::percent(df3$Freq/df3$group)
    }
    else {
      df3$ratio <- scales::percent(df3$ratio1)
    }
    df3$end = cumsum(df3$Freq)
    df3
    df3$start = dplyr::lag(df3$end)
    df3$start[1] = 0
    df3$start1 = df3$start * 2 * pi/total
    df3$end1 = df3$end * 2 * pi/total
    df3$start1 = df3$start1 + start
    df3$end1 = df3$end1 + start
    df3$mid = (df3$start1 + df3$end1)/2
    df3$focus = 0
    if (!is.null(selected)) {
      df3$focus[selected] = explodePos
    }
    else if (!is.null(explode)) {
      selected = c()
      for (i in 1:length(explode)) {
        start = 1 + nrow(a) * (explode[i] - 1)
        selected = c(selected, start:(start + nrow(a) -
                                        1))
      }
      selected
      df3$focus[selected] = explodePos
    }
    df3
    df3$x = 0
    df3$y = 0
    df
    if (!is.null(explode)) {
      explode
      for (i in 1:length(explode)) {
        xpos = df$focus[explode[i]] * sin(df$mid[explode[i]])
        ypos = df$focus[explode[i]] * cos(df$mid[explode[i]])
        df3$x[df3$pie == explode[i]] = xpos
        df3$y[df3$pie == explode[i]] = ypos
      }
    }
    df3$no = 1:nrow(df3)
    df3$label = df3[[donuts]]
    if (showRatioDonut) {
      if (max(nchar(levels(df3$label))) <= 2)
        df3$label = paste0(df3$label, "(", df3$ratio,
                           ")")
      else df3$label = paste0(df3$label, "\n(", df3$ratio,
                              ")")
    }
    df3$label[df3$ratio1 == 0] = ""
    df3$label[df3$ratio1 < showRatioThreshold] = ""
    df3$hjust = ifelse((df3$mid%%(2 * pi)) > pi, 1, 0)
    df3$vjust = ifelse(((df3$mid%%(2 * pi)) < (pi/2)) | (df3$mid%%(2 *
                                                                     pi) > (pi * 3/2)), 0, 1)
    df3$no = factor(df3$no)
    df3
    labelposition
    if (labelposition > 0) {
      df3$radius = r2
      if (explodeDonut)
        df3$radius[df3$focus != 0] = df3$radius[df3$focus !=
                                                  0] + df3$focus[df3$focus != 0]
      df3$segx = df3$radius * sin(df3$mid) + df3$x
      df3$segy = df3$radius * cos(df3$mid) + df3$y
      df3$segxend = (df3$radius + 0.05) * sin(df3$mid) +
        df3$x
      df3$segyend = (df3$radius + 0.05) * cos(df3$mid) +
        df3$y
      if (labelposition == 2)
        df3$radius = (r1 + r2)/2
      df3$labelx = (df3$radius) * sin(df3$mid) + df3$x
      df3$labely = (df3$radius) * cos(df3$mid) + df3$y
    }
    else {
      df3$radius = (r1 + r2)/2
      if (explodeDonut)
        df3$radius[df3$focus != 0] = df3$radius[df3$focus !=
                                                  0] + df3$focus[df3$focus != 0]
      df3$labelx = df3$radius * sin(df3$mid) + df3$x
      df3$labely = df3$radius * cos(df3$mid) + df3$y
    }
    df3$segx[df3$ratio1 == 0] = 0
    df3$segxend[df3$ratio1 == 0] = 0
    df3$segy[df3$ratio1 == 0] = 0
    df3$segyend[df3$ratio1 == 0] = 0
    if (labelposition == 0) {
      df3$segx[df3$ratio1 < showRatioThreshold] = 0
      df3$segxend[df3$ratio1 < showRatioThreshold] = 0
      df3$segy[df3$ratio1 < showRatioThreshold] = 0
      df3$segyend[df3$ratio1 < showRatioThreshold] = 0
    }
    df3
    del = which(df3$Freq == 0)
    del
    if (length(del) > 0)
      subColor <- subColor[-del]
    subColor
  }
  p <- ggplot() + theme_void() + coord_fixed()
  if (is.null(maxx)) {
    r3 = r2 + 0.3
  }
  else {
    r3 = maxx
  }
  p1 <- p + ggforce::geom_arc_bar(aes_string(x0 = "x", y0 = "y",
                                             r0 = as.character(r0), r = as.character(r1), start = "start1",
                                             end = "end1", fill = pies), alpha = pieAlpha, color = color,
                                  data = df) + transparent() + scale_fill_manual(values = mainCol) +
    xlim(r3 * c(-1, 1)) + ylim(r3 * c(-1, 1)) + guides(fill = FALSE)
  if ((labelposition == 1) & (is.null(donuts))) {
    p1 <- p1 + geom_segment(aes_string(x = "segx",
                                       y = "segy", xend = "segxend", yend = "segyend"),
                            data = df) + geom_text(aes_string(x = "segxend",
                                                              y = "segyend", label = "label", hjust = "hjust",
                                                              vjust = "vjust"), size = pieLabelSize, data = df,
                                                   family = family)
  }
  else if ((labelposition == 2) & (is.null(donuts))) {
    p1 <- p1 + geom_segment(aes_string(x = "segx",
                                       y = "segy", xend = "segxend", yend = "segyend"),
                            data = df[df$ratio < labelpositionThreshold, ]) +
      geom_text(aes_string(x = "segxend", y = "segyend",
                           label = "label", hjust = "hjust",
                           vjust = "vjust"), size = pieLabelSize,
                data = df[df$ratio < labelpositionThreshold,
                ], family = family) + geom_text(aes_string(x = "labelx",
                                                           y = "labely", label = "label"), size = pieLabelSize,
                                                data = df[df$ratio >= labelpositionThreshold, ],
                                                family = family)
  }
  else {
    p1 <- p1 + geom_text(aes_string(x = "labelx", y = "labely",
                                    label = "label"), size = pieLabelSize, data = df,
                         family = family)
  }
  if (showPieName)
    p1 <- p1 + annotate("text", x = 0, y = 0, label = pies,
                        size = titlesize, family = family)
  p1 <- p1 + theme(text = element_text(family = family))
  if (!is.null(donuts)) {
    if (explodeDonut) {
      p3 <- p + ggforce::geom_arc_bar(aes_string(x0 = "x",
                                                 y0 = "y", r0 = as.character(r1), r = as.character(r2),
                                                 start = "start1", end = "end1", fill = "no",
                                                 explode = "focus"), alpha = donutAlpha,
                                      color = color, data = df3)
    }
    else {
      p3 <- p + ggforce::geom_arc_bar(aes_string(x0 = "x",
                                                 y0 = "y", r0 = as.character(r1), r = as.character(r2),
                                                 start = "start1", end = "end1", fill = "no"),
                                      alpha = donutAlpha, color = color, data = df3)
    }
    p3 <- p3 + transparent() + scale_fill_manual(values = subColor) +
      xlim(r3 * c(-1, 1)) + ylim(r3 * c(-1, 1)) + guides(fill = FALSE)
    p3
    if (labelposition == 1) {
      p3 <- p3 + geom_segment(aes_string(x = "segx",
                                         y = "segy", xend = "segxend", yend = "segyend"),
                              data = df3) + geom_text(aes_string(x = "segxend",
                                                                 y = "segyend", label = "label", hjust = "hjust",
                                                                 vjust = "vjust"), size = donutLabelSize,
                                                      data = df3, family = family)
    }
    else if (labelposition == 0) {
      p3 <- p3 + geom_text(aes_string(x = "labelx",
                                      y = "labely", label = "label"), size = donutLabelSize,
                           data = df3, family = family)
    }
    else {
      p3 <- p3 + geom_segment(aes_string(x = "segx",
                                         y = "segy", xend = "segxend", yend = "segyend"),
                              data = df3[df3$ratio1 < labelpositionThreshold,
                              ]) + geom_text(aes_string(x = "segxend",
                                                        y = "segyend", label = "label", hjust = "hjust",
                                                        vjust = "vjust"), size = donutLabelSize,
                                             data = df3[df3$ratio1 < labelpositionThreshold,
                                             ], family = family) + geom_text(aes_string(x = "labelx",
                                                                                        y = "labely", label = "label"), size = donutLabelSize,
                                                                             data = df3[df3$ratio1 >= labelpositionThreshold,
                                                                             ], family = family)
    }
    if (!is.null(title))
      p3 <- p3 + annotate("text", x = 0, y = r3,
                          label = title, size = titlesize, family = family)
    else if (showDonutName)
      p3 <- p3 + annotate("text", x = (-1) * r3,
                          y = r3, label = donuts, hjust = 0, size = titlesize,
                          family = family)
    p3 <- p3 + theme(text = element_text(family = family))
    grid::grid.newpage()
    print(p1, vp = grid::viewport(height = 1, width = 1))
    print(p3, vp = grid::viewport(height = 1, width = 1))
  }
  else {
    p1
  }
}


# Pie Donuts --------------------------------------------------------------

meta <- metadata %>%
  dplyr::filter(approach!="/" & sub_type=="trend") %>%
  dplyr::mutate(approach=factor(approach,levels=c("top-down","bottom-up")),
                ccl2=factor(ccl2,levels=c("decreasing trends","it depends trend","mixed trends","increasing trends")))

meta$ccl <- ifelse(meta$ccl2=="it depends trend","factor-\ndependent\ntrends",as.character(meta$ccl2))

meta2 <- metadata %>%
  dplyr::filter(approach!="/" & sub_type=="drivers") %>%
  dplyr::mutate(approach=factor(approach,levels=c("top-down","bottom-up")),
                ccl2=factor(ccl2,levels=c("negative effect of drivers","it depends drivers","no/positive effect of drivers")))
meta2$ccl1 <- ifelse(meta2$ccl2=="it depends drivers","factor-\ndependent\ndrivers",as.character(meta2$ccl2))
meta2$ccl <- ifelse(meta2$ccl1=="no/positive effect of drivers","none/positive\neffect of drivers",as.character(meta2$ccl1))

pie_trends <- PieDonutCustom(meta,aes(pies=approach,donuts=ccl),
                             start=pi/2,maxx=2,
                             r0=0.1,r1=0.8,r2=1.4,showPieName=FALSE,
                             pieLabelSize = 6, donutLabelSize = 4,
                             family = "Montserrat",labelpositionThreshold = 0.03,
                             palette_name = "Dark2"
                             )


pie_drivers <- PieDonutCustom(meta2,aes(pies=approach,donuts=ccl),
                              start=2.5*pi,maxx=2,
                              r0=0.1,r1=0.8,r2=1.4,showPieName=FALSE,
                              pieLabelSize = 6, donutLabelSize = 4,
                              family="Montserrat",
                              showRatioThreshold = 0.01,
                              labelposition = 0,
                              palette_name = "Dark2"
                              )


###############################################################################
###=========================== IMPACT OF DATA ==============================###
###############################################################################


# Create dataframe --------------------------------------------------------
dataf_biodiv <- as.data.frame(table(stats$data,stats$sub_type,stats$taxo,stats$time,stats$ccl2)) %>%
  # select relevant lines
  dplyr::filter(Freq>0 & Var5!="irrelevant" & Var4!="time unknown") %>%
  # rename the columns
  dplyr::rename(data="Var1",
                analysis="Var2",
                group="Var3",
                time="Var4",
                ccl="Var5",
                n="Freq") %>%
  # order factors
  dplyr::mutate(data=factor(data,levels=c("Biotime","LPD","Biotime, LPD","GPDD","Other","Aggregation")),
                group=factor(group,levels=c("4 + groups","2-3 groups","1 group")),
                time=factor(time,levels=c("51 + years","31 - 50 years","- 30 years")),
                ccl=factor(ccl,levels=c("decreasing trends","factor-dependent trend","mixed trends","increasing trends",
                                        "negative effect of drivers","factor-dependent drivers","none/positive effect of drivers"))
  )

# Trends alluvial ---------------------------------------------------------
cols <- hcl.colors(4, "Geyser", rev=T)
(f4_A <- ggplot(data = dataf_biodiv[which(dataf_biodiv$analysis=="trend"),],
               aes(axis1 = data,
                   axis2 = group,
                   axis3 = time,
                   axis4 = ccl,
                   y=n)) +
  stat_alluvium(aes(fill = ccl), alpha=.6, lode.guidance = "backward") +
  geom_stratum(#aes(fill=analysis),
    alpha=.65) +
  scale_fill_manual(values=cols,na.value = "white") +

  # name data
  annotate ("text",x = 1,y = 28,label = "Data",size=4,family="Montserrat",fontface="bold.italic") +
  annotate ("text",x = 1,y = 26,label = "Biotime",size=4,family="Montserrat",fontface="bold") +
  annotate ("text",x = 1,y = 22.5,label = "LPD",size=4,family="Montserrat",fontface="bold") +
  annotate ("text",x = 1,y = 19,label = "Biotime +\nLPD",size=3,family="Montserrat",fontface="bold") +
  annotate ("text",x = 1,y = 17,label = "Other",size=4,family="Montserrat",fontface="bold") +
  annotate ("text",x = 1,y = 9,label = "Aggregation",size=4,family="Montserrat",fontface="bold",angle=90) +

  # name groups
  annotate ("text",x = 2,y = 28,label = "Taxonomic\ngroups",size=4,family="Montserrat",fontface="bold.italic") +
  annotate ("text",x = 2,y = 20,label = "4 + groups",size=4,family="Montserrat",fontface="bold",angle=90) +
  annotate ("text",x = 2,y = 11.5,label = "2-3\ngroups",size=4,family="Montserrat",fontface="bold") +
  annotate ("text",x = 2,y = 5,label = "1 group",size=4,family="Montserrat",fontface="bold",angle=90) +

  # name time
  annotate ("text",x = 3,y = 28,label = "Temporal\nextent",size=4,family="Montserrat",fontface="bold.italic") +
  annotate ("text",x = 3,y = 22.5,label = "50 +\nyears",size=4,family="Montserrat",fontface="bold") +
  annotate ("text",x = 3,y = 14,label = "30 - 49\nyears",size=4,family="Montserrat",fontface="bold") +
  annotate ("text",x = 3,y = 5,label = "- 30\nyears",size=4,family="Montserrat",fontface="bold") +

  # name conclusions
  annotate ("text",x = 4,y = 28,label = "Conclusion",size=4,family="Montserrat",fontface="bold.italic") +
  annotate ("text",x = 4,y = 20,label = "Decreasing trends",size=4,angle=90,family="Montserrat",fontface="bold") +
  annotate ("text",x = 4,y = 11,label = "Factor-\ndependent\ntrends",size=3,family="Montserrat",fontface="bold") +
  annotate ("text",x = 4,y = 5,label = "Mixed\ntrends",size=4,family="Montserrat",fontface="bold") +
  annotate ("text",x = 4,y = 1,label = "Increasing\ntrends",size=3,family="Montserrat",fontface="bold") +

  # visualization parameters
  theme_void() +
  theme(legend.position = "none"))

# Drivers alluvial --------------------------------------------------------
cols <- hcl.colors(3, "Geyser",rev=T)
(f5_A <- ggplot(data = dataf_biodiv[which(dataf_biodiv$analysis=="drivers"),],
               aes(axis1 = data,
                   axis2 = group,
                   axis3 = time,
                   axis4 = ccl,
                   y=n)) +
  stat_alluvium(aes(fill = ccl), alpha=.6, lode.guidance = "backward") +
  geom_stratum(#aes(fill=analysis),
    alpha=.65) +
  scale_fill_manual(values=cols,na.value = "white") +

  # name data
  annotate ("text",x = 1,y = 9.5,label = "Data",size=4,family="Montserrat",fontface="bold.italic") +
  annotate ("text",x = 1,y = 8.5,label = "Biotime",size=4,family="Montserrat",fontface="bold") +
  annotate ("text",x = 1,y = 7,label = "Biotime +\nLPD",size=3,family="Montserrat",fontface="bold") +
  annotate ("text",x = 1,y = 5.5,label = "GPDD",size=4,family="Montserrat",fontface="bold") +
  annotate ("text",x = 1,y = 4.5,label = "Other",size=4,family="Montserrat",fontface="bold") +
  annotate ("text",x = 1,y = 2.25,label = "Aggregation",size=4,family="Montserrat",fontface="bold",angle=90) +

  # name groups
  annotate ("text",x = 2,y = 9.5,label = "Taxonomic\ngroups",size=4,family="Montserrat",fontface="bold.italic") +
  annotate ("text",x = 2,y = 5.5,label = "4 + groups",size=4,family="Montserrat",fontface="bold",angle=90) +
  annotate ("text",x = 2,y = 1.5,label = "2-3\ngroups",size=4,family="Montserrat",fontface="bold") +
  annotate ("text",x = 2,y = 0.5,label = "1\ngroup",size=4,family="Montserrat",fontface="bold") +

  # name time
  annotate ("text",x = 3,y = 9.5,label = "Temporal\nextent",size=4,family="Montserrat",fontface="bold.italic") +
  annotate ("text",x = 3,y = 7,label = "50 +\nyears",size=4,family="Montserrat",fontface="bold") +
  annotate ("text",x = 3,y = 5.5,label = "30 - 49\nyears",size=4,family="Montserrat",fontface="bold") +
  annotate ("text",x = 3,y = 3,label = "- 30\nyears",size=4,family="Montserrat",fontface="bold") +

  # name conclusions
  #geom_segment(aes(x = 4, y = 10, xend = 4.25, yend = 10),size=0.25) +
  annotate ("text",x = 4,y = 9.5,label = "Conclusion",size=4,family="Montserrat",fontface="bold.italic") +
  annotate ("text",x = 4,y = 6,label = "Negative\neffect",size=3,family="Montserrat",fontface="bold") +
  annotate ("text",x = 4,y = 2.5,label = "Factor\ndependent\ndrivers",size=3,family="Montserrat",fontface="bold") +
  annotate ("text",x = 4,y = 1,label = "None/\nPositive\neffect",size=3,family="Montserrat",fontface="bold") +

  # visualization parameters
  theme_void() +
  theme(legend.position = "none"))


###############################################################################
###========================== IMPACT OF METHODS ============================###
###############################################################################


# Trends alluvial ---------------------------------------------------------
cols <- hcl.colors(4, "Geyser",rev = T)
(f4_B <- ggplot(data = methest[which(methest$analysis=="trend"),],
                aes(#axis1 = data,
                  axis2= level,
                  axis3 = method,
                  axis4=ccl,
                  y=n)) +
    stat_alluvium(aes(fill = ccl), alpha=.7, lode.guidance = "backward") +
    geom_stratum(#aes(fill=analysis),
      alpha=.65) +
    scale_fill_manual(values=cols,na.value = "white") +

    # name methods
    annotate ("text",x = 2,y = 33,label = "Method",size=4,family="Montserrat",fontface="bold.italic") +
    annotate ("text",x = 2,y = 26,label = "Linear\nmodels\non\nindividual\ntime\nseries",size=3,family="Montserrat",fontface="bold") +
    annotate ("text",x = 2,y = 13.5,label = "Global\nindicators",size=3,family="Montserrat",fontface="bold") +
    annotate ("text",x = 2,y = 3,label = "Other",size=3,family="Montserrat",fontface="bold") +

    # name level"species","population","species, community","population, species, community"
    annotate ("text",x = 1,y = 33,label = "Ecological\nlevel",size=4,family="Montserrat",fontface="bold.italic") +
    annotate ("text",x = 1,y = 22,label = "Population",size=3,family="Montserrat",fontface="bold") +
    annotate ("text",x = 1,y = 8,label = "Species",size=3,family="Montserrat",fontface="bold") +
    annotate ("text",x = 1,y = 2,label = "Grouped",size=3,family="Montserrat",fontface="bold") +

    # name conclusions
    annotate ("text",x = 3,y = 33,label = "Conclusion",size=4,family="Montserrat",fontface="bold.italic") +
    annotate ("text",x = 3,y = 25,label = "Decreasing trends",size=4,angle=90,family="Montserrat",fontface="bold") +
    annotate ("text",x = 3,y = 15.5,label = "Factor\ndependent\ntrends",size=3,family="Montserrat",fontface="bold") +
    annotate ("text",x = 3,y = 7,label = "Mixed\ntrends",size=4,family="Montserrat",fontface="bold") +
    annotate ("text",x = 3,y = 1,label = "Increasing\ntrends",size=3,family="Montserrat",fontface="bold") +

    # visualization parameters
    theme_void() +
    theme(legend.position = "none"))


# Drivers alluvial --------------------------------------------------------
cols <- hcl.colors(3, "Geyser",rev = T)
(f5_B <- ggplot(data = methest[which(methest$analysis=="drivers"),],
                  aes(#axis1 = data,
                    axis2= level,
                    axis3 = method,
                    axis4=ccl,
                    y=n)) +
    stat_alluvium(aes(fill = ccl), alpha=.7, lode.guidance = "backward") +
    geom_stratum(alpha=.65) +
    scale_fill_manual(values=cols,na.value = "white") +

    # name level"species","population","species, community","population, species, community"
    annotate ("text",x = 1,y = 21,label = "Ecological\nlevel",size=4,family="Montserrat",fontface="bold.italic") +
    annotate ("text",x = 1,y = 17,label = "Population",size=3,family="Montserrat",fontface="bold") +
    annotate ("text",x = 1,y = 10,label = "Species",size=3,family="Montserrat",fontface="bold") +
    annotate ("text",x = 1,y = 3,label = "Grouped",size=3,family="Montserrat",fontface="bold") +

    # name methods
    annotate ("text",x = 2,y = 21,label = "Driver",size=4,family="Montserrat",fontface="bold.italic") +
    annotate ("text",x = 2,y = 15,label = "Climate change",size=4,family="Montserrat",fontface="bold",angle=90) +
    annotate ("text",x = 2,y = 6,label = "Anthropogenic\ndrivers",size=4,family="Montserrat",fontface="bold",angle=90) +
    annotate ("text",x = 2,y = 1,label = "Conservation",size=3,family="Montserrat",fontface="bold") +

    # name conclusions
    #geom_segment(aes(x = 4, y = 10, xend = 4.25, yend = 10),size=0.25) +
    annotate ("text",x = 3,y = 21,label = "Conclusion",size=4,family="Montserrat",fontface="bold.italic") +
    annotate ("text",x = 3,y = 14,label = "Negative\neffect",size=3,family="Montserrat",fontface="bold") +
    annotate ("text",x = 3,y = 5,label = "Factor\ndependent\ndrivers",size=3,family="Montserrat",fontface="bold") +
    annotate ("text",x = 3,y = 1,label = "None/\nPositive\neffect",size=3,family="Montserrat",fontface="bold") +

    # visualization parameters
    theme_void() +
    theme(legend.position = "none"))


###############################################################################
###============================ FINAL FIGURES ==============================###
###############################################################################

(ftrends <- plot_grid(f4_A,f4_B,labels=c("A","B"),ncol=2,nrow=1))
(fdrivers <- plot_grid(f5_A,f5_B,labels=c("A","B"),ncol=2,nrow=1))

###############################################################################
###=========================== SAVING PLOTS ================================###
###############################################################################

## Save to PNG

# Appendix S2
ggsave("publication_year.png", plot = plot_publication_year,
       type = 'cairo',
       width = 8, height = 6, dpi = 150)

# Figure 2
ggsave("data_caracteristics.jpg", plot = fig2,
       type = 'cairo',
       width = 15, height = 8, dpi = 150)

# Figure 4
ggsave("trends.png", plot = ftrends,
       type = 'cairo',
       width = 16, height = 10, dpi = 150)

# Figure 5
ggsave("drivers.png", plot = fdrivers,
       type = 'cairo',
       width = 16, height = 10, dpi = 150)
