"""Actual author checkpoint inference. Standard dictionary plugin, no reference argument."""
import os, sys, types
from pathlib import Path
import numpy as np
import torch
xf = types.ModuleType('xformers')
ops = types.ModuleType('xformers.ops')

def forbidden(*a, **k):
    raise RuntimeError('Unconfigured xformers attention path used')
ops.memory_efficient_attention = forbidden
xf.ops = ops
sys.modules['xformers'] = xf
sys.modules['xformers.ops'] = ops
from omegaconf import OmegaConf
from model.network import initialize_diff_net, initialize_control_net
from model.diffusion import initialize_diff_model
from model.diffusion.gaussian_diffusion import get_named_beta_schedule

class Plugin:

    def __init__(self, source, weights, device='cuda'):
        self.device = device
        torch.set_num_threads(4)
        self.config = OmegaConf.load(str(Path(source) / 'configs/test/test_normal.yaml'))
        self.config.exp.batch_size = 1
        self.config.exp.num_workers = 0
        self.config.net.unet_activation = 'SiLU'
        self.config.diffusion.respacing = 100
        self.model = initialize_diff_net(self.config).to(device).eval()
        self.control = initialize_control_net(self.config).to(device).eval()
        self.model.load_state_dict(torch.load(Path(weights) / 'ToothCraft_normal_diff_branch.pth', map_location='cpu', weights_only=True)['state_dict'], strict=True)
        self.control.load_state_dict(torch.load(Path(weights) / 'ToothCraft_normal_control_branch.pth', map_location='cpu', weights_only=True)['state_dict'], strict=True)
        betas = get_named_beta_schedule(self.config.diffusion.beta_schedule, self.config.diffusion.step, self.config.diffusion.scale_ratio)
        self.diff = initialize_diff_model(betas, self.config, self.config.diffusion.model)
        assert self.diff.num_timesteps == 10, 'Registered 10 transitions required'

    def generate(self, task):
        if set(task) != {'tsdf', 'seed'}:
            raise ValueError('Plugin input schema: no hidden reference allowed')
        x = np.asarray(task['tsdf'], np.float32)
        if x.shape != (64, 64, 64) or not np.isfinite(x).all():
            raise ValueError('Invalid TSDF')
        torch.manual_seed(int(task['seed']))
        torch.cuda.manual_seed_all(int(task['seed']))
        with torch.inference_mode():
            hint = torch.from_numpy(np.clip(x, -1, 1)).to(self.device)[None, None]
            out = self.diff.p_sample_loop(model=self.model, control_model=self.control, shape=[1, 1, 64, 64, 64], device=self.device, progress=True, noise=None, clip_denoised=False, model_kwargs={'hint': hint, 'y': None, 'antag': None, 'noise_save_path': None})
        return dict(status='DESIGN', tsdf=out[0, 0].cpu().numpy(), units='normalized', participant='ToothCraft_normal_author_weights')
