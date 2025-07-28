from django.shortcuts import render
from django.http import JsonResponse
from django.urls import reverse
from django.utils.http import urlsafe_base64_decode
import requests
from pathlib import Path
import json, os, glob
import pandas as pd
import numpy as np
import tempfile
from datetime import datetime
import subprocess
import time
import shutil
from urllib.parse import unquote
from concurrent.futures import ProcessPoolExecutor, as_completed

from apps import getKey, outputPath, outputSamplePath, sampleList, nowTime, multiSampleList
from apps.protein.services import setSample, complexSample, runDownloadPDB, proteinMakeResultFolder, runOpenfold, runAlphafold2, runEsmfold, runBoltzChai, proteinAlignment
# Create your views here.

def proteinInputPageSample(request):
    databases = ['Uniref30_2302', 'PDB70_220313', 'colabfold_envdb_202108']
    samples = setSample()
    # print(samples)
    return render(request, 'protein/input_sample.html', {'databases': databases, 'samples_json': json.dumps(samples)})
        
def proteinProcessSample(request):
    if request.method == "POST":    
        try:
            data = json.loads(request.body)
            sample_name = data.get('sampleName', '')
            results_dir = f"{outputPath()}/protein/_{sample_name}"

            return JsonResponse({"status": "done", "results_dir" : results_dir})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)
        
def proteinOutputPage(request):
    results_dir = request.GET.get('results_dir')
        
    # with open(f"{results_dir}/rot-origin.pdb", 'r') as f1:
    #     origin = f1.read()
        
    pdb_map = {
        "origin_pdb": "rot-origin.pdb",
        "openfold_pdb": "rot-openfold.pdb",
        "af2_pdb": "rot-alphafold.pdb",
        "esmfold_pdb": "rot-esmfold.pdb",
        "boltz_pdb": "rot-boltz.pdb",
        "chai_pdb": "rot-chai.pdb",
    }

    context = {
        key: open(f"{results_dir}/{filename}", 'r').read()
        for key, filename in pdb_map.items()
    }

    return render(request, 'protein/result.html', context)
    
def chaiInputPage(request):
    samples = setSample()
    # print(samples)
    return render(request, 'chai/input.html', {'samples_json': json.dumps(samples)})

def chaiProcess(request):
    if request.method == "POST":    
        try:
            results_dir = f"{outputPath()}/chai"
            
            data = json.loads(request.body)
            sample_name = data.get('sampleName', '')    
            time.sleep(1)  # 5초 동안 멈춤
            return JsonResponse({"status": "done", "results_dir" : results_dir, "sample_name" : sample_name})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)

def chaiOutputPage(request):
    sample_name = request.GET.get('sample_name')
    results_dir = request.GET.get('results_dir')
    
    print(sample_name)
    sample_pdb_path = os.path.join(results_dir, f"{sample_name}.pdb")

    pdb_content = ""
    if os.path.exists(sample_pdb_path):
        with open(sample_pdb_path, 'r') as f:
            pdb_content = f.read()

    return render(request, 'chai/result.html', {
        "result_pdbs": pdb_content
    })

def chaiComplexInputPage(request):
    samples = complexSample()
    # print(samples)
    return render(request, 'chai/complex_input.html', {'samples_json': json.dumps(samples)})

def chaiComplexProcess(request):
    if request.method == "POST":    
        try:
            results_dir = f"{outputPath()}/chai_complex"
            
            data = json.loads(request.body)
            sample_name = data.get('sampleName', '')    
            time.sleep(1.5)  # 5초 동안 멈춤
            return JsonResponse({"status": "done", "results_dir" : results_dir, "sample_name" : sample_name})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)

def chaiComplexOutputPage(request):
    sample_name = request.GET.get('sample_name')
    results_dir = request.GET.get('results_dir')

    # print(sample_name)
    sample_json_path = os.path.join(results_dir, f"{sample_name}.json")

    pdb_content = ""
    if os.path.exists(sample_json_path):
        with open(sample_json_path, 'r') as f:
            pdb_content = json.load(f)

    return render(request, 'chai/complex_result.html', {
        "result_pdbs": pdb_content
    })

# boltz
def boltzInputPage(request):
    samples = setSample()
    # print(samples)
    return render(request, 'boltz/input.html', {'samples_json': json.dumps(samples)})

def boltzProcess(request):
    if request.method == "POST":    
        try:
            results_dir = f"{outputPath()}/boltz"
            
            data = json.loads(request.body)
            sample_name = data.get('sampleName', '')    
            time.sleep(1)  # 5초 동안 멈춤
            return JsonResponse({"status": "done", "results_dir" : results_dir, "sample_name" : sample_name})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)

def boltzOutputPage(request):
    sample_name = request.GET.get('sample_name')
    results_dir = request.GET.get('results_dir')
    
    print(sample_name)
    sample_pdb_path = os.path.join(results_dir, f"{sample_name}.pdb")

    pdb_content = ""
    if os.path.exists(sample_pdb_path):
        with open(sample_pdb_path, 'r') as f:
            pdb_content = f.read()

    return render(request, 'boltz/result.html', {
        "result_pdbs": pdb_content
    })


# boltz complex 
def boltzComplexInputPage(request):
    samples = complexSample()
    # print(samples)
    return render(request, 'boltz/complex_input.html', {'samples_json': json.dumps(samples)})

def boltzComplexProcess(request):
    if request.method == "POST":    
        try:
            results_dir = f"{outputPath()}/boltz_complex"
            
            data = json.loads(request.body)
            sample_name = data.get('sampleName', '')    
            print(sample_name)
            time.sleep(1.5)  # 5초 동안 멈춤
            return JsonResponse({"status": "done", "results_dir" : results_dir, "sample_name" : sample_name})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)

def boltzComplexOutputPage(request):
    sample_name = request.GET.get('sample_name')
    results_dir = request.GET.get('results_dir')
    
    # print(sample_name)
    sample_json_path = os.path.join(results_dir, f"{sample_name}.json")

    pdb_content = ""
    if os.path.exists(sample_json_path):
        with open(sample_json_path, 'r') as f:
            pdb_content = json.load(f)

    return render(request, 'boltz/complex_result.html', {
        "result_pdbs": pdb_content
    })

###################################################################################################
def boltz2OutputListPage(request):
    results_dir = request.GET.get('results_dir')
    source = request.GET.get('source')

    with open(f"{results_dir}/boltz2.json", 'r') as f:
        boltz2 = json.load(f)
    if boltz2 is None:
        return JsonResponse({"error": "세션에 저장된 결과가 없습니다."}, status=404)

    cards = []
    if (source == "single") | (source == "multi"):
        tool_name = "Boltz2"
    else:
        tool_name = "Boltz2 complex"
    
    
    for i in range(len(boltz2)):
        title_pretty = f"Model {i + 1}"

        if source == "single" :
            url = f"{reverse('protein:boltz2OutputPage')}?results_dir={results_dir}&num={i}"
        elif source == "multi":
            url = f"{reverse('protein:boltz2MultiOutputPage')}?results_dir={results_dir}&num={i}"
        elif source == "complex_single":
            url = f"{reverse('protein:boltz2ComplexOutputPage')}?results_dir={results_dir}&num={i}"
        elif source == "complex_multi":
            url = f"{reverse('protein:boltz2MultiComplexOutputPage')}?results_dir={results_dir}&num={i}"
            
        cards.append({
            "title": title_pretty,
            "icon": "img/single_module/boltz.png",
            "color": "#198754",
            "url": url,
            "delay": f"{0.1 + i * 0.1:.1f}s"
        })

    return render(request, 'list.html', {"cards": cards, "tool": tool_name})


def boltz2InputPage(request):
    return render(request, 'boltz2/input.html', {'samples_json': json.dumps(sampleList())})

def boltz2Process(request):
    if request.method == "POST":    
        try:
            sample_dir = f"{outputSamplePath()}/boltz2"
            
            data = json.loads(request.body)
            sample_name = data.get("sampleName") 
            
            results_dir = f"{sample_dir}/{sample_name}"
            
            return JsonResponse({"status": "done", "results_dir" : results_dir})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)
    
def boltz2OutputPage(request):
    results_dir = request.GET.get('results_dir')
    num = int(request.GET.get('num'))

    with open(f"{results_dir}/boltz2.json", 'r') as f:
        fold = json.load(f)
    if fold is None:
        return JsonResponse({"error": "세션에 저장된 결과가 없습니다."}, status=404)

    return render(request, 'boltz2/result.html', {"result_pdb" : fold[num], "num" : num + 1})



def boltz2MultiInputPage(request):
    return render(request, 'boltz2/multimer_input.html', {'samples_json': json.dumps(multiSampleList())})

def boltz2MultiProcess(request):
    if request.method == "POST":    
        try:
            sample_dir = f"{outputSamplePath()}/boltz2_multi"
            
            data = json.loads(request.body)
            sample_name = data.get("sampleName") 
            
            results_dir = f"{sample_dir}/{sample_name}"
            
            return JsonResponse({"status": "done", "results_dir" : results_dir})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)
    
    
def boltz2MultiOutputPage(request):
    results_dir = request.GET.get('results_dir')
    num = int(request.GET.get('num'))

    with open(f"{results_dir}/boltz2.json", 'r') as f:
        fold = json.load(f)
    if fold is None:
        return JsonResponse({"error": "세션에 저장된 결과가 없습니다."}, status=404)

    return render(request, 'boltz2/multimer_result.html', {"result_pdb" : fold[num], "num" : num + 1})
    
# boltz2 complex 
def boltz2ComplexInputPage(request):
    return render(request, 'boltz2/complex_input.html', {'samples_json': json.dumps(sampleList())})

def boltz2ComplexProcess(request):
    if request.method == "POST":    
        try:
            sample_dir = f"{outputSamplePath()}/boltz2_complex"
            
            data = json.loads(request.body)
            sample_name = data.get("sampleName") 
            results_dir = f"{sample_dir}/{sample_name}"
            
            return JsonResponse({"status": "done", "results_dir" : results_dir})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)

def boltz2ComplexOutputPage(request):
    results_dir = request.GET.get('results_dir')
    num = int(request.GET.get('num'))

    with open(f"{results_dir}/boltz2.json", 'r') as f:
        fold = json.load(f)
    if fold is None:
        return JsonResponse({"error": "세션에 저장된 결과가 없습니다."}, status=404)

    return render(request, 'boltz2/complex_result.html', {"result_pdb" : fold[num], "num" : num + 1})


def boltz2MultiComplexInputPage(request):
    return render(request, 'boltz2/complex_multimer_input.html')

def boltz2MultiComplexOutputPage(request):
    pdb_path = request.GET.get('pdb_path')  # 예: /static/output/.../model_0.pdb
    
    # 서버 절대경로로 복원
    full_path = os.path.join("/opt/platform/LNPSolution_AI_Smart_Bench", pdb_path.lstrip("/"))

    pdb_content = ""
    if os.path.exists(full_path):
        with open(full_path, 'r') as f:
            pdb_content = f.read()
    else:
        pdb_content = f"❌ 파일이 존재하지 않습니다: {full_path}"

    return render(request, 'boltz2/complex_multimer_result.html', {
                            "result_pdbs": pdb_content
    })