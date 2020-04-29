# from django.conf.urls import patterns, include, url
from django.contrib import admin
from django.conf import settings
# from django.urls import include, path
from django.urls import path
from django.conf.urls.static import static
# from views import *
from ocean import views
import os

admin.autodiscover()
PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))

urlpatterns = [
                path('', views.main_smiles, name='home'), # if you don't have a MarvinJS licence
                path('admin/', admin.site.urls),
                # url(r'^$', main_marvin, name='home'), # if you have a MarvinJS licence

                # path(r'^media/(?P<path>.*)$', 'django.views.static.serve', {'document_root': '%s/media' % PROJECT_ROOT}),
                path('findNeighbourhoodForCompounds/', views.findNeighbourhoodForCompounds, name='findNeighbourhoodForCompounds'),
                path('ajax_get_cmpd_data/', views.getCmpdsForTarget),
                path('png_for_smiles/', views.png_for_smiles, name='png_for_smiles'),
                path('calc', views.calcOceanStatistics),
                path('createfps', views.createAllFPSforAllMolregnos),
                path('report',views.getOceanReportForSmiles),
                       # Uncomment the admin/doc line below to enable admin documentation:
                       # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),
              ] + static(settings.STATIC_URL,
                         document_root=settings.STATIC_ROOT)

if os.path.exists('ocean/custom_urls.py'):
    import custom_urls
    mcu = {k:v for k,v in custom_urls.__dict__.items() if not k.startswith('__')}
    current_locals = locals()
    for entry,value in mcu.items():
        if not entry in current_locals.keys() or \
                        entry in ['urlpatterns']:
            current_locals.update({entry:value})
            print("load custom urls", entry, value)