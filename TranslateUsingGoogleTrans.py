def translate(text, translateTo):

  # Translate the text provided
  return translator.translate(text, dest=translateTo).text
