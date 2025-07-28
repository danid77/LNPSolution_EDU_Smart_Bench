import tempfile, json, os, subprocess, glob
from datetime import datetime
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import SDWriter

from apps import toolPath

def runMolspark(smiles, model_type, strategy, temperature, num_samples):
    # print(smiles, model_type, strategy, temperature, num_samples)
    # 임시 파일 생성
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as smiles_file:
        smiles_file.write(smiles)  # str 그대로 쓰기
        smiles_path = smiles_file.name  # 파일 경로 저장

    print(f"임시 파일 경로: {smiles_path}")

    # 파일이 잘 생성됐는지 확인
    with open(smiles_path, "r") as f:
        print(f.read())

    _, output_file = tempfile.mkstemp(suffix=".csv")
    
    # 기본 JSON 데이터 구조
    reinvent_option = {
        "run_type": "sampling",
        "use_cuda": True,
        "parameters": {
            "model_file": f"/opt/git-tools/REINVENT4/models/{model_type}.prior",
            "smiles_file": smiles_path,
            "sample_strategy": strategy,
            "output_file": output_file,
            "num_smiles": int(num_samples),
            "unique_molecules": True,
            "randomize_smiles": True
        }
    }
    
    # temperature가 None이 아니면 추가
    if temperature is not None:
        reinvent_option["parameters"]["temperature"] = float(temperature)

    print(reinvent_option)
    # 임시 JSON 파일 생성
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".json") as option_json:
        json.dump(reinvent_option, option_json, indent=4)
        option_json_path = option_json.name

    print(f"임시 JSON 파일 경로: {option_json_path}")
    print(output_file)
    
    reinvent = ["/opt/anaconda3/envs/reinvent4/bin/reinvent", "-f", "json", option_json_path]
    process1 = subprocess.run(reinvent, capture_output=True, text=True)
    print("STDOUT:", process1.stdout)
    print("STDERR:", process1.stderr)
    
    os.remove(smiles_path)
    os.remove(option_json_path)
    # os.remove(output_file)
    
    # csv_list = os.listdir(cal_output_dir)
    # csv_file = [file for file in csv_list if file.find(".csv") !=-1][0]
    # cal_output_file = f"{cal_output_dir}/{csv_file}"
    
    return output_file

def runCal(result_df_csv, results_dir):
    chemical_cal = ["/opt/anaconda3/envs/chemplot/bin/python", f"{toolPath()}/chemical_property_cal/mol_cal_property.py",
                    result_df_csv, results_dir]
    
    process2 = subprocess.run(chemical_cal, capture_output=True, text=True)
    print("STDOUT:", process2.stdout)
    print("STDERR:", process2.stderr)
    
    result_csv_file = glob.glob(f"{results_dir}/*.csv")[0]
    return result_csv_file


def csv_to_dict(_csv_file):
    """ CSV 파일을 Pandas DataFrame으로 읽어 딕셔너리 리스트로 변환 """
    try:
        csv_file = pd.read_csv(_csv_file).astype(object)
        csv_file = csv_file.where(pd.notnull(csv_file), None)  # 🔥 `_csv_file` 대신 `csv_file` 사용!
        return csv_file.to_dict(orient='records')
    except Exception as e:
        print(f"CSV 변환 오류: {e}")
        return None
    
def getCanonicalIsomerSmiles(smiles):
    # # 입체화학 선택 여부 확인
    # while True:
    #     stereo_choice = stereo_value
    #     if stereo_choice in ['y', 'n']:
    #         consider_stereochemistry = True if stereo_choice == 'y' else False
    #         break
    #     else:
    #         print("잘못된 입력입니다. y 또는 n을 입력해주세요.")

    # SMILES로부터 분자 생성
    mol = Chem.MolFromSmiles(smiles)
    
    # 유효한 분자인지 확인
    if mol is None:
        print(f"Invalid SMILES: {smiles}")
        return "Invalid SMILES input."
    
    try:
        isomers = list(AllChem.EnumerateStereoisomers(mol))
        canonical_smiles = [Chem.MolToSmiles(isomer, isomericSmiles=True) for isomer in isomers]
        # print(canonical_smiles)
    except Exception as e:
        print(f"Warning: stereochemistry handling failed. Error: {str(e)}")
        canonical_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
        print("isomer 에러")
    
    return canonical_smiles