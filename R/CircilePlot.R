#' @export

circle.plot<- function(data.raw){
  library(dplyr)
  library(purrr)    # map functions (like lapply)
  library(ggplot2)
  library(lazyeval) # interp function
  library(tidyr)
  library(RColorBrewer)
  library(tibble)
  
  if (any(data.raw$test==0)){
    data.raw$test[data.raw$test==0]<-data.raw$value[data.raw$test==0 & data.raw$variable=="lower"]*0.5
  }
  #data.raw$variable<-factor(data.raw$variable,levels(data.raw$variable)[c(2,1)])
  
  data1 <- tbl_df(data.raw[1:(nrow(data.raw)/2),])
  
  
  # function requires 
  rotate_data <- function(data, col, by_col) {
    lev <- levels(data[,by_col][[1]])
    num <- length(lev)
    
    dir <- rep(seq(((num - 1) * 360 / num), 0, length.out = num))
    
    data$dir_ <- map_dbl(1:nrow(data), function(x) {dir[match(data[x,by_col][[1]], lev)]})
    
    #col_num <- match("mpg", colnames(cars))
    #filter_criteria <- interp(~ which_column == col_num, which_column = as.name(col))
    
    expr <- lazyeval::interp(~x, x = as.name(col))
    data <- mutate_(data, .dots = setNames(list(expr), "plotX"))
    data <- mutate_(data, .dots = setNames(list(expr), "plotY"))
    
    data <- data %>%
      mutate(plotX = round(cos(dir_ * pi / 180) * plotX, 2),
             plotY = round(sin(dir_ * pi / 180) * plotY, 2))
    
    data
  } 
  
  circleFun <- function(center=c(0,0), diameter=1, npoints=100, start=0, end=2, filled=TRUE){
    tt <- seq(start*pi, end*pi, length.out=npoints)
    df <- data.frame(
      x = center[1] + diameter / 2 * cos(tt),
      y = center[2] + diameter / 2 * sin(tt)
    )
    if(filled==TRUE) { #add a point at the center so the whole 'pie slice' is filled
      df <- rbind(df, center)
    }
    return(df)
  }
  
  data1 <- rotate_data(data1, "test", "categories")
  data1<-data1[order(data1$dir_,decreasing = T),]
  
  data_fake <- tbl_df(data.raw)
  data_fake <- rotate_data(data_fake, "value", "categories")
  data_fake<-data_fake[order(data_fake$variable,1/data_fake$dir_,decreasing=F),]
  
  
  line_length <- max(data[,c("value","test")] * 1.1)
  lim <- max(data[,c("value","test")] * 1.1)
  
  rl <- data_frame(dir = unique(data1$dir_), l = rep(line_length, length(unique(data1$dir_)))) %>% 
    mutate(plotX = cos(dir * pi / 180) * (l),
           plotY = sin(dir * pi / 180) * (l))
  rl$xend <- 0
  rl$yend <- 0

  lb <- rl
  lb$label <- levels(data1$categories)
  
  circlegrid <- data_frame(dia = seq(lim / 4, 2 * lim, lim / 4))
  circlegrid <- circlegrid %>% 
    mutate(data = map(dia, function(x) {
      df     <- circleFun(diameter = x, filled = FALSE)
      df$lev <- x
      df
    }))
  
  plotcircles <- bind_rows(circlegrid$data)
  plotcircles$lev <- as.factor(plotcircles$lev)
  
  cl <- data_frame(x = as.numeric(levels(plotcircles$lev)), label = as.character(round(x,1)))
  cl <- cl[cl$x <= max(data$test * 1.1),]

  ggplot() + 
    geom_segment(data = rl, aes(x = plotX, xend = xend, y = plotY, yend = yend), colour = "grey50") +
    geom_path   (data = plotcircles, aes(x = x, y = y, group = lev), colour = "grey50") + 
    geom_text   (data = cl, aes(x = x, y = 1, label = label), colour = "grey40") +
    geom_polygon(data = data_fake, aes(y = plotY, x = plotX), fill = "grey70", colour = 'grey70', size = 1, show.legend = T, alpha = 0.8) +
    geom_path   (data = data1[c(1:nrow(data1),1),], aes(y = plotY, x = plotX), colour = 'steelblue3', size = 1.1) +
    geom_point  (data = data1, aes(y = plotY, x = plotX), stat='identity', colour = 'steelblue4', size = 1.1) +
    geom_text   (data = lb, aes(x = plotX, y = plotY, label = label), colour = "black") +
    ylim(-lim, lim) + xlim(-lim, lim) +
    theme(
      axis.text  = element_blank(), 
      axis.title = element_blank(), 
      line       = element_blank(), 
      rect       = element_blank()
    ) + 
    coord_equal()
  

}