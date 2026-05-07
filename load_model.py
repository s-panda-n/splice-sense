import torch
from transformers import AutoTokenizer, AutoModelForSequenceClassification
from peft import PeftModel

def load_model(model_name):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    if model_name == "nt100m":
        MODEL_PATH = "models/checkpoints/nt-100m-splice"
        tokenizer  = AutoTokenizer.from_pretrained(MODEL_PATH, trust_remote_code=True)
        model      = AutoModelForSequenceClassification.from_pretrained(
            MODEL_PATH, trust_remote_code=True
        )

    elif model_name == "nt500m":
        BASE_MODEL = "InstaDeepAI/nucleotide-transformer-v2-500m-multi-species"
        LORA_PATH  = "models/checkpoints/nt-500m-splice"
        tokenizer  = AutoTokenizer.from_pretrained(LORA_PATH, trust_remote_code=True)
        base       = AutoModelForSequenceClassification.from_pretrained(
            BASE_MODEL, trust_remote_code=True, num_labels=3
        )
        model = PeftModel.from_pretrained(base, LORA_PATH)

    model.eval()
    model.to(device)
    return model, tokenizer, device
