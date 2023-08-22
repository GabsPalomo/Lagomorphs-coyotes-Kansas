library(png)
# coyote.s <- readPNG("sil/Coyote_01.png")
# badger.s <- readPNG("sil/Badger01.png")
# sfox.s <- readPNG("sil/Swiftfox01.png")

get_png <- function(filename) {
  grid::rasterGrob(png::readPNG(filename), interpolate = TRUE)
}

coyote.s <- get_png("sil/Coyote_01.png")
badger.s <- get_png("sil/Badger01.png")
sfox.s  <- get_png("sil/Swiftfox01.png")
cotton.s <- get_png("sil/Cottontail01.png")
jack.s <- get_png("sil/Black_tailed_jackrabbit02.png")
np.s <- get_png("sil/no_predators.png")