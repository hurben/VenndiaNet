MLV
======

### Authors

Benjamin Hur, Dongwon Kang, Sangseon Lee, Sun Kim



### Structures

#### ```/var/www/htdocs/MLV_v2```
- MLV : main djangp app
- work : all works are stored in this folder

#### ```/var/www/html/MLV_static```
- all static files are served here. (bootstrap, jquery, fontawesome, img/js/css ...)

### How to run

python manage.py runserver 0.0.0.0:8082

### If static file needs to be updated ...

python manage.py collectstatic
