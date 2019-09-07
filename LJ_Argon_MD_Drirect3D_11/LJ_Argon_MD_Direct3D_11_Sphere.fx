//--------------------------------------------------------------------------------------
// File: Tutorial09.fx
//
// Copyright (c) Microsoft Corporation. All rights reserved.
//--------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------
// Constant Buffer Variables
//--------------------------------------------------------------------------------------
Texture2D txDiffuse : register( t0 );

cbuffer cbNeverChanges : register( b0 )
{
    float3 vLightDir;
};

cbuffer cbChangesEveryFrame : register( b1 )
{
    matrix WorldViewProj;
    matrix World;
    float4 vMeshColor;
};

struct VS_INPUT
{
    float3 Pos          : POSITION;         // position
    float3 Norm         : NORMAL;           // normal
};

struct PS_INPUT
{
    float4 Pos : SV_POSITION;
    float4 Diffuse : COLOR0;
};


//--------------------------------------------------------------------------------------
// Vertex Shader
//--------------------------------------------------------------------------------------
PS_INPUT VS( VS_INPUT input )
{
    PS_INPUT output = (PS_INPUT)0;
    output.Pos = mul( float4(input.Pos,1), WorldViewProj );
    float3 vNormalWorldSpace = normalize( mul( input.Norm, (float3x3)World ) );

    float fLightAmbient = 0.45f;
    float fLighting = saturate( dot( vNormalWorldSpace, vLightDir ) ) + fLightAmbient;
    float fThreshold = 0.75f;
    output.Diffuse.rgb = fLighting > fThreshold ? fThreshold : fLighting;
    output.Diffuse.a = 1.0f; 

    return output;
}


//--------------------------------------------------------------------------------------
// Pixel Shader
//--------------------------------------------------------------------------------------
float4 PS( PS_INPUT input) : SV_Target
{
    //calculate lighting assuming light color is <1,1,1,1>
    return input.Diffuse * vMeshColor;
}
