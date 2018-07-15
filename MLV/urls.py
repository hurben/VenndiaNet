"""MLV URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.10/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.conf.urls import url, include
    2. Add a URL to urlpatterns:  url(r'^blog/', include('blog.urls'))
"""
from django.conf.urls import url
from django.contrib import admin
from . import views

urlpatterns = [
    url(r'^admin/', admin.site.urls),
    url(r'^$', views.main, name='main'),
    url(r'^work/generate/$', views.generate),           # parameter: id
    url(r'^work/(?P<id>\w+)/delete/$', views.delete),   # only works in admin mode
    url(r'^work/(?P<id>\w+)/s1/$', views.step_upload, name='step_upload'),
    url(r'^work/(?P<id>\w+)/s2/$', views.step_select, name='step_select'),
    url(r'^work/(?P<id>\w+)/s3/$', views.step_params, name='step_params'),
    url(r'^work/(?P<id>\w+)/result/$', views.result, name='result')
]
