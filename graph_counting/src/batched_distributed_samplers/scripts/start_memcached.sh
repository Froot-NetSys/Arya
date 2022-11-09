pkill memcached ; ufw allow 11211 ; memcached -d -m 100000 -p 11211 -u root -M -t 4
