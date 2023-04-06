
#Learn in which diretory youa re workin
getwd()

#Check history of codes
history()

#Adjust your directory
setwd("C:/Users/sedat/Downloads/Bismillah-to-Single_Cell")

#install package rwhatsapp
install.packages("rwhatsapp")

#this is the program retrieved from page https://cran.r-project.org/web/packages/rwhatsapp/vignettes/Text_Analysis_using_WhatsApp_data.html

history <- system.file("extdata", "sample.txt", package = "rwhatsapp")

library("rwhatsapp")
library("dplyr")
#assign all ur whatsapp history to item chat
chat <- rwa_read("C:/Users/sedat/Downloads/WhatsApp Chat - Practice w_ Games&Repeats (1)/_chat.txt") %>% 
  filter(!is.na(author)) # remove messages without author
chat

#MESSAGES PER DAY
library("ggplot2"); theme_set(theme_minimal())
library("lubridate")
library(scales)
chat %>%
  mutate(day = date(time)) %>%
  filter(day >= ymd("2023-01-01") & day <= ymd("2023-01-28")) %>%
  count(day) %>%
  ggplot(aes(x = day, y = n)) +
  geom_bar(stat = "identity") +
  ylab("") + xlab("") +
  ggtitle("Messages per day")

#NUMBER OF MESSAGES -ranking of members with bar graph

chat %>%
  mutate(day = date(time)) %>%
  count(author) %>%
  ggplot(aes(x = reorder(author, n), y = n)) +
  geom_bar(stat = "identity") +
  ylab("") + xlab("") +
  coord_flip() + theme(axis.text.y = element_text(size = 6)) +
  ggtitle("Number of messages")

#FAVOURITE EMOJIS
library("ggimage")
emoji_data <- rwhatsapp::emojis %>% # data built into package
  mutate(hex_runes1 = gsub("\\s.*", "", hex_runes)) %>% # ignore combined emojis
  mutate(emoji_url = paste0("https://abs.twimg.com/emoji/v2/72x72/", 
                            tolower(hex_runes1), ".png"))
(THERE IS A PROBLEM IN URL)

chat %>%
  unnest(emoji) %>%
  count(author, emoji, sort = TRUE) %>%
  group_by(author) %>%
  top_n(n = 6, n) %>%
  left_join(emoji_data, by = "emoji") %>% 
  ggplot(aes(x = reorder(emoji, n), y = n, fill = author)) +
  geom_col(show.legend = FALSE) +
  ylab("") +
  xlab("") +
  coord_flip() +
  geom_image(aes(y = n + 20, image = emoji_url)) +
  facet_wrap(~author, ncol = 2, scales = "free_y") +
  ggtitle("Most often used emojis") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

#FAVOURITE WORDS