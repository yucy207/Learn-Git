
library(plotrix)

flower_plot <- function(sample, value, start, a, b, 
	ellipse_col = rgb(135, 206, 235, 150, max = 255), 
	circle_col = rgb(0, 162, 214, max = 255),
	circle_text_cex = 1.5,
	circle_text = NULL,
	clockwise = FALSE,
	alpha     = 0.8,
	...
	) {
	par( bty = "n", ann = F, xaxt = "n", yaxt = "n", mar = c(1,1,1,1))
	plot(c(0,10),c(0,10),type="n")
	n   <- length(sample)
	deg <- 360 / n


	ellipse_col <- unlist(lapply(ellipse_col, function(t){ a <- col2rgb(t); return(rgb(a[1,1], a[2,1], a[3,1], alpha * 255, max = 255))}))
	if (clockwise) {
		sample <- c(sample[1], rev(tail(sample,-1)))
		value  <- c(value[1], rev(tail(value,-1)))
	}
	data <- data.frame(t(sapply(1:n, function(t){
		ellipse_x      <- 5 + cos((start + deg * (t - 1)) * pi / 180)
		ellipse_y      <- 5 + sin((start + deg * (t - 1)) * pi / 180)
		ellipse_rotate <- deg * (t - 1)

		unique_text_x  <- 5 + (b + 0.5) * cos((start + deg * (t - 1)) * pi / 180)
		unique_text_y  <- 5 + (b + 0.5) * sin((start + deg * (t - 1)) * pi / 180)


		return(c(ellipse_x, ellipse_y, ellipse_rotate, 
			     unique_text_x, unique_text_y))
	})))

	colnames(data) <- c("ellipse_x", "ellipse_y", "ellipse_rotate",
	                    "unique_text_x", "unique_text_y")

	draw.ellipse(x = data$ellipse_x, y = data$ellipse_y, col = ellipse_col, a = a, b = b, angle = data$ellipse_rotate)
	text(x = data$unique_text_x, y = data$unique_text_y, value)

	res <- lapply(1:n, function(t){ 
		sample_text_x      <- 5 + (b + 1.3) * cos((start + deg * (t - 1)) * pi / 180)
		sample_text_y      <- 5 + (b + 1.3) * sin((start + deg * (t - 1)) * pi / 180)
		sample_text_rotate <- deg * (t - 1) + start
		adj                <- 0
		if (deg * (t - 1) < 180 && deg * (t - 1) > 0) {
			sample_text_rotate <- deg * (t - 1) - start
			adj                <- 1
		}
		text( x = sample_text_x, y = sample_text_y, sample[[t]], adj = adj, cex = circle_text_cex, srt = sample_text_rotate)
	})

	draw.circle(x = 5, y = 5, r = 1.3 , col = circle_col, border = circle_col)
	text(x=5, y = 5, circle_text)
}
