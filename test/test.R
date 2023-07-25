getwd()
setwd("/home/data/t050446/Non-tumor research/DXMarkers")
# Load the devtools package
library(usethis)
library(utils)
library(devtools)
# install.packages(c("devtools", "roxygen2", "usethis", "testthat"))
#detach(package:DXMarkers)
rm(list = c("annotate_markers"))
# 在R控制台中运行以下命令
library(roxygen2)
roxygen2::roxygenize()



# Check the package
#check()
#devtools::document()

# 假设你的R包项目目录是"DXMarkers"，请将路径替换为你的实际路径
getwd()
setwd("/home/data/t050446/Non-tumor research/DXMarkers")
# 使用devtools包进行打包
devtools::build("/home/data/t050446/Non-tumor research/DXMarkers")


#######生成并查看说明书#####
# 使用usethis包添加一个名为"my-vignette"的vignette
usethis::use_vignette("DXMarkers")

# Load devtools package
library(devtools)

# Build your vignettes
build_vignettes()

# 查看使用说明书
browseVignettes("DXMarkers")
