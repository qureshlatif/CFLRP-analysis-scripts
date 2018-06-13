library(BCRDataAPI)
library(dplyr)

setwd("/home/RMBO.LOCAL/quresh.latif/CFLRP")

BCRDataAPI::reset_api()
BCRDataAPI::set_api_server('192.168.137.180')
BCRDataAPI::add_columns(c('TransectNum|str',
                          'Point|int',
                          'centerPointEasting|int',
                          'centerPointNorthing|int',
                          'centerPointZone|int')
)
BCRDataAPI::filter_on('Stratum in CO-CFLRP-CF,CO-BCR16-RC,CO-BCR16-PC,CO-BCR16-VO,CO-BCR15-PO')
grab <- BCRDataAPI::get_data()
keep <- grab %>% unique
write.csv(keep, "Bird_survey_point_coords.csv", row.names = F)
