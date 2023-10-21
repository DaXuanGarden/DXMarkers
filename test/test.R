getwd()
setwd("/home/data/t050446/DX_Package/DXMarkers")
# Load the devtools package
library(usethis)
library(utils)
library(devtools)
# install.packages(c("devtools", "roxygen2", "usethis", "testthat"))
#detach(package:DXMarkers)
# rm(list = c("annotate_markers"))
# 在R控制台中运行以下命令
library(roxygen2)
roxygen2::roxygenize()


devtools::document()
# Check the package
# check()


# 假设你的R包项目目录是"DXMarkers"，请将路径替换为你的实际路径
getwd()
setwd("/home/data/t050446/DX_Package/DXMarkers")
# 使用devtools包进行打包
devtools::build("/home/data/t050446/DX_Package/DXMarkers")
library(devtools)
devtools::install_local("/home/data/t050446/DX_Package/DXMarkers_1.0.tar.gz", dependencies = TRUE, upgrade = FALSE)

library(DXMarkers)
setwd("/home/data/t050446/01 Single Cell Project/MR&ScRNA")

DXMarkers_heatmap(scRNA, DXMarkers_result, "00DXMarkers.png")

DXMarkers_dotplots(scRNA, "SCINA_markers.csv", "DXMarkers")

DXMarkers_dotplots(scRNA, "~/SCINA_markers.csv", "~/DXMarkers")



#######生成并查看说明书#####
# 使用usethis包添加一个名为"my-vignette"的vignette
usethis::use_vignette("DXMarkers")

# Load devtools package
library(devtools)

# Build your vignettes
build_vignettes()

# 查看使用说明书
browseVignettes("DXMarkers")
######图片配色选择######
if(!require(paletteer))install.packages("paletteer")
if(!require(scico))install.packages('scico')
if(!require(nord))install.packages('nord')
library(paletteer)
paletteer_c("scico::berlin", n = 50)
paletteer_d("RColorBrewer::Paired")
paletteer_dynamic("cartography::green.pal", 5)
paletteer_d("Polychrome::kelly", n = 22)
paletteer_d("pals::kelly")
paletteer_d("palettesForR::Browns")

#remotes::install_github("mtennekes/cols4all")
#install.packages("colorblindcheck")
library(cols4all)
#交互面板
c4a_gui()
#可以通过函数提取配色(色板名称+所需颜色数量)
mycolor <-c4a("muted",9)
######SCP配色######
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
mycolors <- getPalette(60)
mycolors
