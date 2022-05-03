#' Get iso3 codes of recipient countries in COVAX
#'
#' \url{https://www.who.int/initiatives/act-accelerator/covax}
#'
#' @return Vector of iso3 codes.
#' @export
get_covax_iso3c <- function(){
  c(
    "AFG",  "AGO",  "BDI",  "BEN",  "BFA",  "BGD",  "BOL",  "BTN",
    "CAF",  "CIV",  "CMR",  "COD",  "COG",  "COM",  "CPV",
    "DJI",  "DMA",  "DZA",  "EGY",  "ERI",  "ETH",  "FJI",
    "FSM",  "GHA",  "GIN",  "GMB",  "GNB",  "GRD",  "GUY",
    "HND",  "HTI",  "IDN",  "IND",  "KEN",  "KGZ",  "KHM",  "KIR",
    "LAO",  "LBR",  "LCA",  "LKA",  "LSO",  "MAR",  "MDA",  "MDG",
    "MDV",  "MHL",  "MLI",  "MMR",  "MNG",  "MOZ",  "MRT",  "MWI",
    "NER",  "NGA",  "NIC",  "NPL",  "PAK",  "PHL",  "PNG",  "PRK",
    "PSE",  "RWA",  "SDN",  "SEN",  "SLB",  "SLE",  "SLV",
    "SOM",  "SSD",  "STP",  "SWZ",  "SYR",  "TCD",  "TGO",  "TJK",
    "TLS",  "TON",  "TUN",  "TUV",  "TZA",  "UGA",  "UKR",  "UZB",
    "VCT",  "VNM",  "VUT",  "WSM",  "XKX",  "YEM",  "ZMB",  "ZWE"
  )
}
#' Get iso3 codes of recipient countries in GAVI
#'
#' @source \url{https://www.gavi.org/}
#'
#' @return Vector of iso3 codes.
#' @export
get_gavi_iso3c <- function(){
  c(
    "AFG","BGD","BEN","BFA","BDI","KHM","CMR","CAF","TCD","COM","COD","COG",
    "CIV","DJI","ERI","ETH","GMB","GHA","GIN","GNB","HTI","KEN","PRK","KGZ",
    "LAO","LSO","LBR","MDG","MWI","MLI","MRT","MOZ","MMR","NPL","NIC","NER",
    "NGA","PAK","PNG","RWA","STP","SEN","SLE","SLB","SOM","SSD","SDN","SYR",
    "TJK","TZA","TGO","UGA","UZB","YEM","ZMB","ZWE"
  )
}
#' Get iso3 codes of recipient countries in TheGlobalFund
#'
#' @source \url{https://www.theglobalfund.org/en/}
#'
#' @return Vector of iso3 codes.
#' @export
get_global_fund_iso3c <- function(){
  c(
    "AFG", "DZA", "AGO", "ARM", "AZE", "BGD", "BLR", "BLZ", "BEN", "BTN", "BOL",
    "BWA", "BFA", "BDI", "CPV", "KHM", "CAF", "TCD", "COL", "COG", "COD", "CRI",
    "CIV", "CUB", "DJI", "DOM", "ECU", "EGY", "SLV", "SWZ", "ETH", "GAB", "GMB",
    "GEO", "GHA", "GTM", "GIN", "GNB", "GUY", "IDN", "IRN", "KAZ", "KEN", "KGZ",
    "LAO", "LSO", "LBR", "MDG", "MWI", "MLI", "MRT", "MDA", "MNG", "MAR", "MOZ",
    "NAM", "NPL", "NER", "NGA", "PAK", "PAN", "PRY", "PER", "PHL", "ROU", "RWA",
    "SEN", "SLE", "SLB", "SOM", "ZAF", "SSD", "LKA", "SDN", "TJK", "TZA", "THA",
    "TLS", "TGO", "TUN", "TKM", "UGA", "UKR", "UZB", "VEN", "VNM", "ZMB", "ZWE"
  )
}
#' Get the World Bank Income Categories for a given set of iso3 codes
#'
#' Note that some countries that are currently unclassified have been given either
#' their previous classification (i.e Venezuela) or the classification of the country
#' that they are a territory of (i.e. Hong Kong).
#'
#' @source \url{https://datahelpdesk.worldbank.org/knowledgebase/articles/906519-world-bank-country-and-lending-groups}
#' @param iso3cs Vector of iso 3 codes.
#' @return Ordered factor of classifications codes.
#' @export
get_income_group <- function(iso3cs){
  #list categories
  LIC <- c(
    "AFG", "GNB", "SOM", "BFA", "PRK", "SSD", "BDI", "LBR", "SDN", "CAF", "MDG",
    "SYR", "TCD", "MWI", "TGO", "COD", "MLI", "UGA", "ERI", "MOZ", "YEM", "ETH",
    "NER", "GMB", "RWA", "GIN", "SLE"
  )
  LMIC <- c(
    "AGO", "HND", "PHL", "DZA", "IND", "WSM", "BGD", "IDN", "STP", "BLZ", "IRN",
    "SEN", "BEN", "KEN", "SLB", "BTN", "KIR", "LKA", "BOL", "KGZ", "TZA", "CPV",
    "LAO", "TJK", "KHM", "LSO", "TLS", "CMR", "MRT", "TUN", "COM", "FSM", "UKR",
    "COG", "MNG", "UZB", "CIV", "MAR", "VUT", "DJI", "MMR", "VNM", "EGY", "NPL",
    "PSE", "SLV", "NIC", "ZMB", "SWZ", "NGA", "ZWE", "GHA", "PAK", "HTI", "PNG"
  )
  UMIC <- c(
    "ALB", "GAB", "NAM", "ASM", "GEO", "MKD", "ARG", "GRD", "PAN", "ARM", "GTM",
    "PRY", "AZE", "GUY", "PER", "BLR", "IRQ", "ROU", "BIH", "JAM", "RUS", "BWA",
    "JOR", "SRB", "BRA", "KAZ", "ZAF", "BGR", "LCA", "CHN", "LBN", "VCT", "COL",
    "LBY", "SUR", "CRI", "MYS", "THA", "CUB", "MDV", "TON", "DMA", "MHL", "TUR",
    "DOM", "MUS", "TKM", "GNQ", "MEX", "TUV", "ECU", "MDA", "FJI", "MNE", "VEN"
  )
  HIC <- c(
    "AND", "GRC", "POL", "ATG", "GRL", "PRT", "ABW", "GUM", "PRI", "AUS", "HKG",
    "QAT", "AUT", "HUN", "SMR", "BHS", "ISL", "SAU", "BHR", "IRL", "SYC", "BRB",
    "IMN", "SGP", "BEL", "ISR", "SXM", "BMU", "ITA", "SVK", "VGB", "JPN", "SVN",
    "BRN", "KOR", "ESP", "CAN", "KWT", "KNA", "CYM", "LVA", "LIE", "SWE", "CHL",
    "LTU", "CHE", "HRV", "LUX", "TWN", "CUW", "MAC", "TTO", "CYP", "MLT", "TCA",
    "CZE", "MCO", "ARE", "DNK", "NRU", "GBR", "EST", "NLD", "USA", "FRO", "NCL",
    "URY", "FIN", "NZL", "VIR", "FRA", "MNP", "PYF", "NOR", "DEU", "OMN", "GIB",
    "PLW", "GUF"
  )
  #now create ordered factor as output
  output <- factor(dplyr::case_when(
    iso3cs %in% LIC ~ 4,
    iso3cs %in% LMIC ~ 3,
    iso3cs %in% UMIC ~ 2,
    iso3cs %in% HIC ~ 1,
    TRUE ~ as.numeric(NA)
  ),
         levels = seq(4),
         labels = c(
           "HIC",
           "UMIC",
           "LMIC",
           "LIC"
         ), ordered = T)
  if(any(is.na(output))){
    warning(paste0("Unable to get assign the following iso3cs an income group: ",
                   paste0(iso3cs[is.na(output)], collapse = ", ")))
  }
  return(output)
}
#' Get the WHO Regions for a given set of iso3 codes
#'
#' Note that some countries that are currently unclassified have been given the
#' classification of the country that they are a territory of (i.e. Hong Kong).
#'
#' @source \url{https://www.who.int/countries}
#' @param iso3cs Vector of iso 3 codes.
#' @return Vector of WHO regions.
#' @export
get_WHO_region <- function(iso3cs){
  AFR <- c(
    "DZA", "AGO", "BEN", "BWA", "BFA", "BDI", "CPV", "CMR", "CAF", "TCD", "COM",
    "COG", "CIV", "COD", "GNQ", "ERI", "SWZ", "ETH", "GAB", "GMB", "GHA", "GIN",
    "GNB", "KEN", "LSO", "LBR", "MDG", "MWI", "MLI", "MRT", "MUS", "MOZ", "NAM",
    "NER", "NGA", "RWA", "STP", "SEN", "SYC", "SLE", "ZAF", "SSD", "TGO", "UGA",
    "TZA", "ZMB", "ZWE"
  )
  AMR <- c(
    "ATG", "ARG", "BHS", "BRB", "BLZ", "BOL", "BRA", "CAN", "CHL", "COL", "CRI",
    "CUB", "DMA", "DOM", "ECU", "SLV", "GRD", "GTM", "GUY", "HTI", "HND", "JAM",
    "MEX", "NIC", "PAN", "PRY", "PER", "KNA", "LCA", "VCT", "SUR", "TTO", "USA",
    "URY", "VEN", "GUF", "ABW", "CUW"
  )
  SEAR <- c(
    "BGD", "BTN", "PRK", "IND", "IDN", "MDV", "MMR", "NPL", "LKA", "THA", "TLS"
  )
  EUR <- c(
    "ALB", "AND", "ARM", "AUT", "AZE", "BLR", "BEL", "BIH", "BGR", "HRV", "CYP",
    "CZE", "DNK", "EST", "FIN", "FRA", "GEO", "DEU", "GRC", "HUN", "ISL", "IRL",
    "ISR", "ITA", "KAZ", "KGZ", "LVA", "LTU", "LUX", "MLT", "MCO", "MNE", "NLD",
    "MKD", "NOR", "POL", "PRT", "MDA", "ROU", "RUS", "SMR", "SRB", "SVK", "SVN",
    "ESP", "SWE", "CHE", "TJK", "TUR", "TKM", "UKR", "GBR", "UZB"
  )
  EMR <- c(
    "AFG", "BHR", "DJI", "EGY", "IRN", "IRQ", "JOR", "KWT", "LBN", "LBY", "MAR",
    "OMN", "PAK", "QAT", "SAU", "SOM", "SDN", "SYR", "TUN", "ARE", "YEM", "PSE"
  )
  WPR <- c(
    "AUS", "BRN", "KHM", "CHN", "COK", "FJI", "JPN", "KIR", "LAO", "MYS", "MHL",
    "FSM", "MNG", "NRU", "NZL", "NIU", "PLW", "PNG", "PHL", "KOR", "WSM", "SGP",
    "SLB", "TON", "TUV", "VUT", "VNM", "HKG", "TWN", "MAC", "PYF", "NCL"
  )
  #create our output
  output <- dplyr::case_when(
    iso3cs %in% AFR ~ "AFR",
    iso3cs %in% AMR ~ "AMR",
    iso3cs %in% SEAR ~ "SEAR",
    iso3cs %in% EUR ~ "EUR",
    iso3cs %in% EMR ~ "EMR",
    iso3cs %in% WPR ~ "WPR",
    TRUE ~ as.character(NA)
  )
  if(any(is.na(output))){
    warning(paste0("Unable to get assign the following iso3cs a WHO region: ",
                   paste0(iso3cs[is.na(output)], collapse = ", ")))
  }
  return(output)
}
