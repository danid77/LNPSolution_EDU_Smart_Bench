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

from apps import outputPath, nowTime
from apps.protein.services import setSample, complexSample, runDownloadPDB, proteinMakeResultFolder, runOpenfold, runAlphafold2, runEsmfold, runBoltzChai, proteinAlignment
# Create your views here.

def proteinInputPage(request):
    databases = ['Uniref30_2302', 'PDB70_220313', 'colabfold_envdb_202108']
    samples = setSample()
    # print(samples)
    return render(request, 'protein/input.html', {'databases': databases, 'samples_json': json.dumps(samples)})

def proteinInputPageSample(request):
    databases = ['Uniref30_2302', 'PDB70_220313', 'colabfold_envdb_202108']
    samples = setSample()
    # print(samples)
    return render(request, 'protein/input_sample.html', {'databases': databases, 'samples_json': json.dumps(samples)})

# def proteinProcess(request):
#     if request.method == "POST":    
#         try:
#             current_time = datetime.now().strftime("%Y%m%d_%Hh%Mm%Ss")
#             results_dir = f"{outputPath()}/protein/{current_time}"
            
#             data = json.loads(request.body)
#             sample_name = data.get('sampleName', '')
#             protein = data.get('aminoSeq', '')
#             msa_db = data.get('msaDatabases', '')
#             print("Input protein sequence: ", protein)
            
#             os.makedirs(results_dir, exist_ok=True)
#             openfold_result_dir, alphafold_result_dir, esmfold_result_dir, boltz_result_dir, chai_result_dir = proteinMakeResultFolder(results_dir)
#             # print(opnefold_result)
            
#             # runOpenfold(openfold_result_dir, msa_db, protein)
#             # runAlphafold2(alphafold_result_dir, protein)
#             # runEsmfold(esmfold_result_dir, protein)
#             # runBoltzChai(results_dir, boltz_result_dir, chai_result_dir, protein)
            
#             print("Origin protein PDB download.....")
#             runDownloadPDB(results_dir, sample_name)
            
#             # 병렬 실행
#             with ProcessPoolExecutor(max_workers=4) as executor:
#                 futures = {
#                     executor.submit(runOpenfold, openfold_result_dir, msa_db, protein): "Openfold",
#                     executor.submit(runAlphafold2, alphafold_result_dir, protein, sample_name): "Alphafold2",
#                     executor.submit(runEsmfold, esmfold_result_dir, protein): "Esmfold",
#                     executor.submit(runBoltzChai, results_dir, boltz_result_dir, chai_result_dir, protein): "Boltz&Chai"
#                 }

#                 for future in as_completed(futures):
#                     task_name = futures[future]
#                     try:
#                         result = future.result()
#                         print(f"{task_name} 완료: {result}")
#                     except Exception as e:
#                         print(f"{task_name} 실패: {e}")
#                         raise  # 예외를 다시 발생시켜서 전체 중단
            
#             Structure_alignment = proteinAlignment(results_dir)
#             print(Structure_alignment)

#             return JsonResponse({"status": "done", "results_dir" : results_dir})  # ← 클라이언트에 저장 완료만 알림
#         except Exception as e:
#             print("에러 발생:", e)
#             return JsonResponse({"error": str(e)}, status=500)
        
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


def boltz2InputPage(request):
    samples = setSample()
    # print(samples)
    return render(request, 'boltz2/input.html', {'samples_json': json.dumps(samples)})

def boltz2MultiInputPage(request):
    return render(request, 'boltz2/multimer_input.html')

def boltz2Process(request):
    if request.method == "POST":    
        try:
            results_dir = f"{outputPath()}/boltz2/{nowTime()}"
            
            data = json.loads(request.body)
            sample_name = data.get("sampleName")
            sequence = data.get("sequence")
            recycling_steps = int(data.get("recycling_steps"))
            sampling_steps = int(data.get("sampling_steps"))
            num_samples = int(data.get("num_samples"))   
            
            print(sample_name, sequence, recycling_steps, sampling_steps, num_samples)
            
            os.makedirs(results_dir, exist_ok=True)
            # runBoltz2(results_dir, sequence, recycling_steps, sampling_steps, num_samples, None)
            
            return JsonResponse({"status": "done", "results_dir" : results_dir})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)

def boltz2OutputListPage(request):
    results_dir = request.GET.get('results_dir')
    source = request.GET.get('source')

    pdb_paths = glob.glob(f"{results_dir}/*/predictions/*/*.pdb")
    pdb_paths.sort()

    if not pdb_paths:
        return render(request, 'boltz2/list.html', {
            "cards": [],
            "pdb_paths": [],
            "error_message": "❌ 생성된 모델이 없습니다. 입력값을 확인하거나 다시 시도해 주세요."
        })

    cards = []
    for i, pdb_path in enumerate(pdb_paths):
        file_name = os.path.basename(pdb_path)
        title = file_name.replace(".pdb", "")
        title_pretty = f"Model {i + 1}"
        rel_pdb_path = pdb_path.replace("/opt/platform/LNPSolution_AI_Smart_Bench", "")

        if source == "single" :
            url = f"{reverse('protein:boltz2OutputPage')}?pdb_path={rel_pdb_path}"
        elif source == "multi":
            url = f"{reverse('protein:boltz2MultiOutputPage')}?pdb_path={rel_pdb_path}"
        elif source == "complex_single":
            url = f"{reverse('protein:boltz2ComplexOutputPage')}?pdb_path={rel_pdb_path}"
        elif source == "complex_multi":
            url = f"{reverse('protein:boltz2MultiComplexOutputPage')}?pdb_path={rel_pdb_path}"
            
        cards.append({
            "title": title_pretty,
            "icon": "img/single_module/boltz.png",
            "color": "#198754",
            "url": url,
            "delay": f"{0.1 + i * 0.1:.1f}s"
        })

    return render(request, 'boltz2/list.html', {
        "cards": cards,
        "pdb_paths": pdb_paths
    })
    
def boltz2OutputPage(request):
    pdb_path = request.GET.get('pdb_path')  # 예: /static/output/.../model_0.pdb

    # 서버 절대경로로 복원
    full_path = os.path.join("/opt/platform/LNPSolution_AI_Smart_Bench", pdb_path.lstrip("/"))

    pdb_content = ""
    if os.path.exists(full_path):
        with open(full_path, 'r') as f:
            pdb_content = f.read()
    else:
        pdb_content = f"❌ 파일이 존재하지 않습니다: {full_path}"

    return render(request, 'boltz2/result.html', {
        "result_pdbs": pdb_content
    })
    
def boltz2MultiOutputPage(request):
    pdb_path = request.GET.get('pdb_path')  # 예: /static/output/.../model_0.pdb

    # 서버 절대경로로 복원
    full_path = os.path.join("/opt/platform/LNPSolution_AI_Smart_Bench", pdb_path.lstrip("/"))

    pdb_content = ""
    if os.path.exists(full_path):
        with open(full_path, 'r') as f:
            pdb_content = f.read()
    else:
        pdb_content = f"❌ 파일이 존재하지 않습니다: {full_path}"

    return render(request, 'boltz2/multimer_result.html', {
        "result_pdbs": pdb_content
    })
    
# boltz2 complex 
def boltz2ComplexInputPage(request):
    return render(request, 'boltz2/complex_input.html')

def boltz2MultiComplexInputPage(request):
    return render(request, 'boltz2/complex_multimer_input.html')

def boltz2ComplexProcess(request):
    if request.method == "POST":    
        try:
            results_dir = f"{outputPath()}/boltz2_complex/{nowTime()}"
            
            data = json.loads(request.body)
            # sample_name = data.get("sampleName")
            sequence = data.get("sequence")
            ligands = data.get("ligands")
            recycling_steps = int(data.get("recycling_steps"))
            sampling_steps = int(data.get("sampling_steps"))
            num_samples = int(data.get("num_samples"))   
            
            print("sequence : ", sequence)
            print("ligands : ", ligands)
            print(recycling_steps, sampling_steps, num_samples)
            
            os.makedirs(results_dir, exist_ok=True)
            runBoltz2(results_dir, sequence, recycling_steps, sampling_steps, num_samples, ligands)
            
            return JsonResponse({"status": "done", "results_dir" : results_dir})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)

def boltz2ComplexOutputPage(request):
    pdb_path = request.GET.get('pdb_path')  # 예: /static/output/.../model_0.pdb
    
    # 서버 절대경로로 복원
    full_path = os.path.join("/opt/platform/LNPSolution_AI_Smart_Bench", pdb_path.lstrip("/"))

    pdb_content = ""
    if os.path.exists(full_path):
        with open(full_path, 'r') as f:
            pdb_content = f.read()
    else:
        pdb_content = f"❌ 파일이 존재하지 않습니다: {full_path}"

    return render(request, 'boltz2/complex_result.html', {
                            "result_pdbs": pdb_content
    })

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