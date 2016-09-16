from django.conf.urls import patterns, include, url
from views import *
import os
from django.contrib import admin

admin.autodiscover()
PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))

urlpatterns = patterns('',
                       url(r'^admin/', include(admin.site.urls)),
                       # url(r'^$', main_marvin, name='home'), # if you have a MarvinJS licence
                       url(r'^$', main_smiles, name='home'), # if you don't have a MarvinJS licence
                       url(r'^media/(?P<path>.*)$', 'django.views.static.serve', {'document_root': '%s/media' % PROJECT_ROOT}),
                       url(r'^findNeighbourhoodForCompounds/$', findNeighbourhoodForCompounds, name='findNeighbourhoodForCompounds'),
                       url(r'^ajax_get_cmpd_data/$', getCmpdsForTarget),
                       url(r'^png_for_smiles/$', png_for_smiles, name='png_for_smiles'),
                       url(r'^calc/$', calcOceanStatistics),
                       url(r'^createfps/$',createAllFPSforAllMolregnos),
                       # Uncomment the admin/doc line below to enable admin documentation:
                       # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),
                       )

if os.path.exists('ocean/custom_urls.py'):
    import custom_urls
    mcu = {k:v for k,v in custom_urls.__dict__.items() if not k.startswith('__')}
    current_locals = locals()
    for entry,value in mcu.items():
        if not entry in current_locals.keys() or \
                        entry in ['urlpatterns']:
            current_locals.update({entry:value})
            print "load custom urls",entry,value