#!/bin/sh
cd /home/fix/Desktop
while true
do
  /home/fix/Desktop/ngrok  tcp 22 -region=eu > /dev/null &
  sleep 5
  curl  http://localhost:4040/api/tunnels > tunnels.json
  echo "###########################url_address"
  curl  http://localhost:4040/api/tunnels
  echo "#######################################"
  cat testing tunnels.json > sent
  /usr/sbin/sendmail -t < sent
  sleep 12h
  killall ngrok
  echo "&&&&&&&&&  hi xiao, ngrok rebooted!!!! &&&&&&&&&&"
done


