#' Data set of Ozone Data
#'
#' 112 daily measurements of meteorological variables (windspeed, temperature, rainfall, cloudiness)
#' and ozone concentration recorded in Rennes (France) during the summer 2001.
#' In this study, an individual is a day. Thirteen variables were measured: 11 numerical variables
#'  and 2 categorical variables with 2 or 4l evels. There is no missing values.
#' In more detail, the variables are:
#'\describe{
#'  \item{maxO3}{maximum of daily ozone concentration measured in gr.m-3}
#'  \item{T9, T12, T15}{daily temperatures measured in degree Celsius at 9, 12 and 15h (numerical variables called temperature variables hereafter)}
#'  \item{Ne9, Ne12, Ne15}{cloudiness measured at 9, 12 and 15h (numerical variablescalled cloudiness variables hereafter)}
#'  \item{Vx9, Vx12, Vx15}{wind speed (E-W component) measured at 9, 12 and 15h (numerical variables called wind variables hereafter)}
#'  \item{maxO3v}{maximum concentration of ozone measured the day before}
#'  \item{vent}{wind direction measured at 12h (categorical variable with four categories)}
#'  \item{pluie}{occurrence or not of rainfall (categorical variable with two categories)}
#' }
#' The initial objective is to explain the maximum of daily ozone concentration (the response variable is \emph{thusmaxO3}) by the
#' 10 numerical variables available (T9, T12, T15, Ne9, Ne12, Ne15, Vx9, Vx12, Vx15 and maxO3v).
#' The twocategorical variables \emph{vent} and \emph{pluie} are not used in this study.
#'
#' @usage
#' data(ozone)
#'
#' @references{
#'   \insertRef{cornillon2012r}{outlierSIR}
#' }
#'
#' @importFrom Rdpack reprompt
#'
"ozone"
