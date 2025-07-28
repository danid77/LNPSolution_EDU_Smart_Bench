from django.urls import path

from . import views

app_name = 'nim'

urlpatterns = [
    path('nimOpenfoldInputPage/', views.nimOpenfoldInputPageSample, name='nimOpenfoldInputPageSample'),
    path('nimOpenfoldProcess/', views.nimOpenfoldProcessSample, name='nimOpenfoldProcessSample'),
    path('nimOpenfoldList/', views.nimOpenfoldList, name='nimOpenfoldList'),
    path('nimOpenfoldOutputPage/', views.nimOpenfoldOutputPage, name='nimOpenfoldOutputPage'),
    
    path('nimGenmolInputPage/', views.nimGenmolInputPageSample, name='nimGenmolInputPageSample'),
    path('nimGenmolProcess/', views.nimGenmolProcessSample, name='nimGenmolProcessSample'),
    path('nimGenmolOutputPage/', views.nimGenmolOutputPage, name='nimGenmolOutputPage'),
    
    path('nimMolminInputPage/', views.nimMolminInputPageSample, name='nimMolminInputPageSample'),
    path('nimMolminProcess/', views.nimMolminProcessSample, name='nimMolminProcessSample'),
    path('nimMolminOutputPage/', views.nimMolminOutputPage, name='nimMolminOutputPage'),
    
    path('nimAlphafoldInputPage/', views.nimAlphafoldInputPage, name='nimAlphafoldInputPage'),
    path('nimAlphafoldProcess/', views.nimAlphafoldProcess, name='nimAlphafoldProcess'),
    path('nimAlphafoldList/', views.nimAlphafoldList, name='nimAlphafoldList'),
    path('nimAlphafoldOutputPage/', views.nimAlphafoldOutputPage, name='nimAlphafoldOutputPage'),
    
    path('nimAlphafoldMultiInputPage/', views.nimAlphafoldMultiInputPage, name='nimAlphafoldMultiInputPage'),
    path('nimAlphafoldMultiProcess/', views.nimAlphafoldMultiProcess, name='nimAlphafoldMultiProcess'),
    path('nimAlphafoldMultiListPage/', views.nimAlphafoldMultiListPage, name='nimAlphafoldMultiListPage'),
    path('nimAlphafoldMultiOutputPage/', views.nimAlphafoldMultiOutputPage, name='nimAlphafoldMultiOutputPage'),
    
    path('nimEsmfoldInputPage/', views.nimEsmfoldInputPage, name='nimEsmfoldInputPage'),
    path('nimEsmfoldProcess/', views.nimEsmfoldProcess, name='nimEsmfoldProcess'),
    path('nimEsmfoldList/', views.nimEsmfoldList, name='nimEsmfoldList'),
    path('nimEsmOutputPage/', views.nimEsmOutputPage, name='nimEsmOutputPage'),
    
    path('nimRfdiffusionInputPage/', views.nimRfdiffusionInputPage, name='nimRfdiffusionInputPage'),
    path('nimRfdiffusionProcess/', views.nimRfdiffusionProcess, name='nimRfdiffusionProcess'),
    path('nimRfdiffusionOutputPage/', views.nimRfdiffusionOutputPage, name='nimRfdiffusionOutputPage'),
    
    path('nimDiffdockInputPage/', views.nimDiffdockInputPage, name='nimDiffdockInputPage'),
    path('nimDiffdockProcess/', views.nimDiffdockProcess, name='nimDiffdockProcess'),
    path('nimDiffdockList/', views.nimDiffdockList, name='nimDiffdockList'),
    path('nimDiffdockOutputPage/', views.nimDiffdockOutputPage, name='nimDiffdockOutputPage'),
    
    path('nimProteinmpnnInputPage/', views.nimProteinmpnnInputPage, name='nimProteinmpnnInputPage'),
    path('nimProteinmpnnProcess/', views.nimProteinmpnnProcess, name='nimProteinmpnnProcess'),
    path('nimProteinmpnnList/', views.nimProteinmpnnList, name='nimProteinmpnnList'),
    path('nimProteinmpnnOutputPage/', views.nimProteinmpnnOutputPage, name='nimProteinmpnnOutputPage'),
]