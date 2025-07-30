from django.shortcuts import render
from django.http import JsonResponse
from django.http import HttpResponse
import json, os, tempfile, subprocess
from datetime import datetime
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# Create your views here.
from apps import outputPath, nowTime, outputSamplePath
from apps.generatemol import services as sv
# Create your views here.

def molsparkInputPage(request):
    return render(request, 'molspark/input.html')

def molsparkProcess(request):
    try:
        sample_dir = f"{outputSamplePath()}/molspark"
        
        # JSON 데이터 파싱
        data = json.loads(request.body)
        molSeq = data.get("molSeq")
        strategy = data.get("strategy")
        
        print(molSeq, strategy)
        
        results_dir = f"{sample_dir}/{molSeq}_{strategy}"
        
        return JsonResponse({"status": "success", "results_dir": results_dir})
    except Exception as e:
        return JsonResponse({"status": "error", "message": str(e)})
    
# def molsparkOutputPage(request):
#     return render(request, 'molspark/output.html')

def molsparkOutputPage(request):
    results_dir = request.GET.get('results_dir')
    if not results_dir or not os.path.isdir(results_dir):
        return JsonResponse({'error': 'Invalid results_dir'}, status=400)
    
    # genmol_df = pd.read_csv(f"{results_dir}/genmol_result.csv")
    result_csv_file = pd.read_csv(f"{results_dir}/Rdkit_calculate_result.csv").drop(columns=[f"fp_{i+1}" for i in range(2048)])
    
    # 현재 컬럼 목록
    current_columns = result_csv_file.columns.tolist()

    # 앞으로 보낼 컬럼
    front_cols = ['Status', 'SMILES']

    # 나머지 컬럼에서 front_cols 제거
    remaining_cols = [col for col in current_columns if col not in front_cols]

    # 새 순서로 재정렬
    new_column_order = front_cols + remaining_cols

    # 컬럼 순서 재배치
    result_csv_file = result_csv_file[new_column_order]
    result_csv_file = result_csv_file.reset_index().rename(columns={"index" : "Num"})
    
    return render(request, 'molspark/output.html', {
                                                      "smiles" : result_csv_file['SMILES'].tolist(), 
                                                      "cal_data": result_csv_file.replace({np.nan: "Null"}).to_dict(orient="records")
                                                      })
    
def chemDrawInputPage(request):
    return render(request, 'chemdraw/input.html', 
                  {"ketcher_path" : "tools/streamlit_ketcher/streamlit_ketcher/frontend/index.html"})
    
def chemDrawProcess(request):
    try:
        results_dir = f"{outputPath()}/chemdraw/{nowTime()}"
        if not os.path.exists(results_dir):
            os.mkdir(results_dir)
        
        # JSON 데이터 파싱
        data = json.loads(request.body)
        smiles_input = data.get("smilesInput")
        
        # print(smiles_input)
        
        canonical_smiles = sv.getCanonicalIsomerSmiles(smiles_input)
        print("canonical_smiles : ", canonical_smiles)
        pd.DataFrame({"SMILES" : canonical_smiles}).to_csv(f"{results_dir}/isomer_SMILES.csv", index=False)
        
        return JsonResponse({"status": "success", "results_dir": results_dir})
    except Exception as e:
        return JsonResponse({"status": "error", "message": str(e)})
    
def chemDrawOutputPage(request):
    results_dir = request.GET.get('results_dir')
    if not results_dir or not os.path.isdir(results_dir):
        return JsonResponse({'error': 'Invalid results_dir'}, status=400)
    
    # genmol_df = pd.read_csv(f"{results_dir}/genmol_result.csv")
    result_csv_file = pd.read_csv(f"{results_dir}/isomer_SMILES.csv")

    return render(request, 'chemdraw/result.html',{
        "smiles" : result_csv_file['SMILES'].tolist()
    })
    
def chemDrawSmilesTo3d(request):
    try:
        body = json.loads(request.body)
        smiles = body.get("smiles")
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(mol)
        mol_block = Chem.MolToMolBlock(mol)
        return JsonResponse({"status": "success", "mol_block": mol_block})
    except Exception as e:
        return JsonResponse({"status": "error", "message": str(e)}, status=500)
    
def pocket2MolInputPage(request):
    samples = [
        {"protein" : "TNIK", "pdbid": "2X7F", "x_coor" : "27.31", "y_coor" : "-5.47", "z_coor" : "57.1"},
        {"protein" : "HER2", "pdbid": "8VB5", "x_coor" : "-1.13", "y_coor" : "6.73", "z_coor" : "-18.11"},
        {"protein" : "BCL2", "pdbid": "2W3L", "x_coor" : "39.81", "y_coor" : "26.94", "z_coor" : "-12.41"},
    ]
    return render(request, 'pocket2mol/input.html', {
        'samples_json': json.dumps(samples)
    })
    
def pocket2MolProcess(request):
    try:
        sample_dir = f"{outputSamplePath()}/pocket2mol"
        
        # JSON 데이터 파싱
        data = json.loads(request.body)
        sample_name = data.get("pdbId")
        
        # print(sample_name)
        
        results_dir = f"{sample_dir}/{sample_name}"
        
        return JsonResponse({"status": "success", "results_dir": results_dir})
    except Exception as e:
        return JsonResponse({"status": "error", "message": str(e)})

def pocket2MolOutputPage(request):
    results_dir = request.GET.get('results_dir')
    if not results_dir or not os.path.isdir(results_dir):
        return JsonResponse({'error': 'Invalid results_dir'}, status=400)

    result_csv_file = pd.read_csv(f"{results_dir}/Rdkit_result.csv")
    
    return render(request, 'pocket2mol/result.html', {"smiles" : result_csv_file['SMILES'].tolist(), 
                                                      "cal_data" : result_csv_file.replace({np.nan: "Null"}).to_dict(orient="records")})